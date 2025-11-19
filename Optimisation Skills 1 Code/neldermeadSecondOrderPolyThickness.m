% ---------------------------------------------------------
% Second-Order (Quadratic) Thickness Optimisation
% ---------------------------------------------------------

clc; clear all; close all;

engFuncs = makeEngFuncs;

% --- Fixed Parameters ---
chord = 0.15;
dist_from_neutral_axis = 0.12 * chord;
N = 100;
targetStress = -30e6; % Pa
penaltyFactor = 1e7;

% --- VALID Initial Guess ---
% a = start thickness, b = linear term, c = quadratic term
X0 = [0.02, -0.005, 0];  % Safe feasible start

% --- BOX CONSTRAINTS ---
a_min = 0.001;  a_max = 0.05;
b_min = -0.05;  b_max = 0.05;
c_min = -0.05;  c_max = 0.05;

lb = [a_min, b_min, c_min];
ub = [a_max, b_max, c_max];

% --- Objective Function Handle ---
objective = @(X) objectiveFunctionQuad(X, chord, dist_from_neutral_axis, N, targetStress, penaltyFactor, engFuncs);

% --- Optimisation Options ---
options = optimset('Display','iter', 'TolX',1e-4, 'TolFun',1e-6);

% --- Run Bounded Nelder-Mead ---
[X_opt, fval] = fminsearchbnd(objective, X0, lb, ub, options);

% --- Extract Optimised Parameters ---
a_opt = X_opt(1);
b_opt = X_opt(2);
c_opt = X_opt(3);

% --- Recompute final design ---
[maxStress_opt, volume_opt] = evaluateDesignQuad(a_opt, b_opt, c_opt, chord, dist_from_neutral_axis, N, engFuncs);

fprintf('\n*** Optimization Complete ***\n');
fprintf('Optimal a (root thickness): %.4f m\n', a_opt);
fprintf('Optimal b (linear term): %.4f\n', b_opt);
fprintf('Optimal c (quadratic term): %.4f\n', c_opt);
fprintf('Minimum Volume: %.6e m^3\n', volume_opt);
fprintf('Maximum Bending Stress: %.3f MPa\n', maxStress_opt / 1e6);
fprintf('Target Stress (magnitude): %.0f MPa\n', abs(targetStress) / 1e6);

%% --------------------- PLOTS ---------------------
figure;
theme(gcf, "light");

subplot(2,1,1);
dists = linspace(0,1,N);
thicknesses = a_opt + b_opt*dists + c_opt*dists.^2;
plot(dists, thicknesses, 'LineWidth', 2);
xlabel('Normalised Radial position, r/R');
ylabel('Thickness, t');
title("Optimised Quadratic Thickness Distribution");

subplot(2,1,2);
I_array_opt = engFuncs.secondMomentAreaArraySkinThicknessQuadratic(N, chord, a_opt, b_opt, c_opt);
[~, M, ~, ~, ~] = beamBending(N, I_array_opt);
bendingStresses_opt = engFuncs.bendingStress(M, dist_from_neutral_axis, I_array_opt);

hold on;
plot(dists, abs(bendingStresses_opt)/1e6, 'b-', 'LineWidth', 2);
line([0,1], [abs(targetStress)/1e6, abs(targetStress)/1e6], 'Color','r','LineStyle','--');
xlabel('Normalised Radial position, r/R');
ylabel('Bending Stress, MPa');
legend('Computed |Stress|', 'Target Stress', 'Location', 'best');
grid on;
hold off;

%% ---------------- Objective Function ----------------
function [f, maxStress, volume] = objectiveFunctionQuad(X, chord, dist_from_neutral_axis, N, targetStress, penaltyFactor, engFuncs)
    a = X(1); b = X(2); c = X(3);

    % --- Check thickness non-negative at critical points ---
    % Evaluate at x=0, x=1, and extremum if c ~= 0
    x_ext = -b/(2*c);
    t0 = a; t1 = a + b + c;
    t_ext = Inf;
    if c ~= 0 && x_ext >= 0 && x_ext <= 1
        t_ext = a + b*x_ext + c*x_ext^2;
    end
    t_min = min([t0, t1, t_ext]);
    if t_min <= 0
        f = Inf; maxStress = Inf; volume = Inf;
        return;
    end

    % --- Evaluate design ---
    [maxStress, volume] = evaluateDesignQuad(a, b, c, chord, dist_from_neutral_axis, N, engFuncs);

    % --- Stress penalty ---
    targetMag = abs(targetStress);
    maxStressMag = abs(maxStress);
    if targetMag == 0
        stressRatio = 0;
    else
        stressRatio = (maxStressMag - targetMag)/targetMag;
    end
    penalty = max(0, penaltyFactor*stressRatio^2);

    % --- Objective function ---
    f = volume + penalty;

    fprintf('a=%.4f, b=%.4f, c=%.4f -> Volume: %.6e, MaxStress: %.3f MPa, Obj: %.6e\n', ...
        a, b, c, volume, maxStressMag/1e6, f);
end

%% ---------------- Evaluate Design ----------------
function [maxStress, volume] = evaluateDesignQuad(a, b, c, chord, dist_from_neutral_axis, N, engFuncs)
    % --- Evaluate thickness at critical points ---
    x = linspace(0,1,N);
    t_array = a + b*x + c*x.^2;
    
% Ensure thickness never negative
if any(t_array <= 0)
    maxStress = Inf; volume = Inf;
    return;
end

% Ensure monotonic decreasing
x = linspace(0,1,N);
t_derivative = b + 2*c*x;
if any(t_derivative > 0)
    maxStress = Inf; volume = Inf;
    return;
end


    % --- Second moment of area array ---
    I_array = engFuncs.secondMomentAreaArraySkinThicknessQuadratic(N, chord, a, b, c);
    if any(I_array <= 0)
        maxStress = Inf; volume = Inf;
        return;
    end

    % --- Compute volume ---
    volume = engFuncs.findVolumeSkinThicknessQuadratic(chord, a, b, c);

    % --- Beam bending and stresses ---
    [~, M, ~, ~, ~] = beamBending(N, I_array);
    bendingStresses = engFuncs.bendingStress(M, dist_from_neutral_axis, I_array);

    maxStress = max(abs(bendingStresses));
end
