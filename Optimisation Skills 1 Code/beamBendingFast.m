% beam bending fast
function [delta, M, V, L_vec, P]  = beamBendingFast(N, I)

    %% 1. Define Geometry, Material, and Load (Unchanged)
    R = 1; 
    r = linspace(0, R, N); 
    dr = r(2) - r(1); 
    dr4 = dr^4;
    
    E = 100e9; 
    I = I(:); 
    EI = E * I;
    
    L_func = @(r_norm) 50 .* r_norm .* (1 - (1 ./ exp(-50 .* (r_norm - R))));
    L_vec = L_func(r); 
    F = L_vec * dr;
    P = zeros(1, N); 

    %% 2. Construct the Stiffness Matrix K (Unchanged)
    K = zeros(N, N);
    
    % Interior Nodes (i = 3 to N-2)
    for i = 3:N-2 
        K(i, i-2) = EI(i) * 1;    
        K(i, i-1) = EI(i) * (-4); 
        K(i, i)   = EI(i) * 6;    
        K(i, i+1) = EI(i) * (-4); 
        K(i, i+2) = EI(i) * 1;    
    end
    
    % Boundary Node N-1 
    K(N-1, N-3) = EI(N-1) * 1;
    K(N-1, N-2) = EI(N-1) * (-4);
    K(N-1, N-1) = EI(N-1) * 6;
    K(N-1, N)   = EI(N-1) * (-4); 
    
    % Boundary Node N (Tip)
    K(N, N-2) = EI(N) * 1;
    K(N, N-1) = EI(N) * (-2); 
    K(N, N)   = EI(N) * 1;
    
    K = K / dr4;
    
    %% 3. Apply Cantilever Boundary Conditions and Solve for Delta (Unchanged)
    
    % Solve the reduced system: K_mod * delta_mod = F_mod
    K_mod = K(3:N, 3:N);
    F_mod = F(3:N)'; 
    
    % Solve the linear system
    delta_mod = K_mod \ F_mod; 
    
    % Reconstruct the full deflection vector: delta(1)=0, delta(2)=0
    delta = [0; 0; delta_mod]; 
    delta = delta'; 

    %% 4. Calculate Internal Forces (M and V) via Integration (FIX)
    
    % The internal forces must satisfy the free tip conditions: M(R)=0, V(R)=0.
    % We integrate backwards from the tip (R) to the root (0).

    % Initialize V and M vectors
    V = zeros(1, N);
    M = zeros(1, N);
    
    % Tip boundary conditions: V(N) = 0 and M(N) = 0
    V(N) = 0; 
    M(N) = 0;

    % Loop backwards from N-1 down to 1
    for i = N-1:-1:1
        % Shear Force V(i) = V(i+1) + Integral(L_vec) between i and i+1
        % Using the Trapezoidal Rule: Integral(L) ≈ (L_i + L_{i+1}) / 2 * dr
        V(i) = V(i+1) + (L_vec(i) + L_vec(i+1)) * 0.5 * dr;

        % Bending Moment M(i) = M(i+1) + Integral(V) between i and i+1
        % Using the Trapezoidal Rule: Integral(V) ≈ (V_i + V_{i+1}) / 2 * dr
        % Note: We must use the just-calculated V(i) in the trapezoid sum.
        M(i) = M(i+1) + (V(i) + V(i+1)) * 0.5 * dr; 
    end
    
    % V and M are now calculated using a stable integration process.
    % M(1) will now automatically be the maximum reaction moment.

end