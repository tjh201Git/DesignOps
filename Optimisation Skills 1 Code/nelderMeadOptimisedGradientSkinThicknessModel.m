% Nelder mead optimised gradient skin thickness

% --- Nelder-Mead Optimization Script ---

clc, clear all, close all

engFuncs = makeEngFuncs;

% --- Fixed Parameters (from your original script) ---
chord = 0.15;
dist_from_neutral_axis = 0.12 * chord;
N = 100;
targetStress = 30e6; % Pa (30 MPa)
penaltyFactor = 1e7; % Hefty penalty for violating the stress constraint

% --- Initial Guess for [startThickness, gradient] ---
% Based on your original script's values:
X0 = [0.05, -0.3];

% --- Run Nelder-Mead Optimization ---
% Define the anonymous objective function handle
objective = @(X) objectiveFunction(X, chord, dist_from_neutral_axis, N, targetStress, penaltyFactor, engFuncs);

% Use fminsearch (Nelder-Mead)
options = optimset('Display','iter', 'TolX', 1e-4, 'TolFun', 1e-6); % Set optimization options
[X_opt, fval] = fminsearch(objective, X0, options);

% --- Output Results ---
startThickness_opt = X_opt(1);
gradient_opt = X_opt(2);

% Re-evaluate the optimal design to get final stress and volume
[maxStress_opt, volume_opt] = evaluateDesign(startThickness_opt, gradient_opt, chord, dist_from_neutral_axis, N, engFuncs);

fprintf('\n*** Optimization Complete ***\n');
fprintf('Optimal Start Thickness: %.4f m\n', startThickness_opt);
fprintf('Optimal Gradient: %.4f\n', gradient_opt);
fprintf('Minimum Volume: %.6e m^3\n', volume_opt);
fprintf('Maximum Bending Stress: %.3f MPa\n', maxStress_opt / 1e6);
fprintf('Target Stress: %.0f MPa\n', targetStress / 1e6);



%% plots

% 1. Re-evaluate the stress and I array using the optimal parameters
I_array_opt = engFuncs.secondMomentAreaArraySkinThicknessGradient(N, chord, startThickness_opt, gradient_opt);
[~, M, ~, ~, ~] = beamBending(N, I_array_opt); % Requires beamBending function
bendingStresses_opt = engFuncs.bendingStress(M, dist_from_neutral_axis, I_array_opt);

% 2. Create distance array
dists = linspace(0, chord, N);

% 3. Plotting Bending Stress
figure;
theme(gcf, "light");
hold on;
plot(dists, abs(bendingStresses_opt)/1e6, 'b-', 'LineWidth', 2);
line([0, chord], [targetStress/1e6, targetStress/1e6], 'Color', 'r', 'LineStyle', '--');

xlabel('Distance along Chord ($x$), m', 'Interpreter', 'latex');
ylabel('Absolute Bending Stress ($|\sigma|$), MPa', 'Interpreter', 'latex');
title('Bending Stress Distribution along the Aerofoil Chord');
legend('Computed Stress', 'Target Stress (30 MPa)', 'Location', 'best');
grid on;
hold off;

% --- Helper Function to compute thickness along the chord ---
% This re-uses the logic from getThicknessAtDistanceFromLinearLine

function t_at_x = computeThicknessProfile(x_vector, startThickness, gradient)
    t_at_x = zeros(size(x_vector));
    for i = 1:length(x_vector)
        distance = x_vector(i);
        tempThickness = startThickness + gradient * distance;
        positiveThickness = max(0, tempThickness);
        % Ensures thickness doesn't exceed the start thickness
        t_at_x(i) = min(positiveThickness, startThickness); 
    end
end


% 4. Plotting 3D Aerofoil and Thickness 
L_z = chord * 2; % Extruded length for visualization

% 4a. Chord coordinates
x_chord = linspace(0, chord, N);
% NACA0012 profile height (y coordinate)
y_profile = engFuncs.findAerofoilPositiveHeight(x_chord); 

% 4b. Variable thickness profile
thickness_profile = computeThicknessProfile(x_chord, startThickness_opt, gradient_opt);

% 4c. Setup Z-plane coordinates for the 3D plot
Z = [0, L_z]; % Two points in the Z-direction for visualization

figure;
theme(gcf, "light");
hold on;
view(3); % Set to 3D view
axis equal;

% --- Plot the Outer Aerofoil Surface ---
% Top Surface
[X_top, Z_top] = meshgrid(x_chord, Z);
Y_top = repmat(y_profile, 2, 1);
surf(X_top, Z_top, Y_top, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Bottom Surface
Y_bottom = repmat(-y_profile, 2, 1);
surf(X_top, Z_top, Y_bottom, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% --- Plot the Internal Cavity (Showing Thickness) ---
% The skin is defined by the profile y_p and thickness t_p.
% The internal surface y_i is y_p - t_p/2 (assuming t is perpendicular to the surface, simplified here)
% For visualization, let's just plot the thickness at the max height (dist_from_neutral_axis)

% Plot the inner cavity surface (simplified by offsetting the outer surface by the thickness)
% Top Inner Surface
Y_inner_top = repmat(y_profile - thickness_profile/2, 2, 1); 
surf(X_top, Z_top, Y_inner_top, 'FaceColor', [0.1 0.1 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.8); 

% Bottom Inner Surface
Y_inner_bottom = repmat(-y_profile + thickness_profile/2, 2, 1); 
surf(X_top, Z_top, Y_inner_bottom, 'FaceColor', [0.1 0.1 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.8);

% Labeling and Aesthetics
xlabel('Chord ($x$) Distance (m)', 'Interpreter', 'latex');
ylabel('Extrusion ($z$) Distance (m)', 'Interpreter', 'latex');
zlabel('Height ($y$) (m)', 'Interpreter', 'latex');
title('3D Aerofoil Profile with Variable Skin Thickness');
colormap(flipud(gray)); % Use a nice colormap

% Add markers for clarity
plot3(x_chord(1), L_z/2, y_profile(1), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
text(x_chord(1), L_z/2, y_profile(1)*1.2, ['t_0 = ', num2str(startThickness_opt, '%.4f')], 'Color', 'r');

hold off;


% 1. Chord coordinates
x_chord = linspace(0, chord, N);

% 2. NACA0012 profile height (y coordinate)
y_profile = engFuncs.findAerofoilPositiveHeight(x_chord); 

% 3. Variable thickness profile
thickness_profile = computeThicknessProfile(x_chord, startThickness_opt, gradient_opt);

% 4. Define the inner and outer boundaries for the skin visualization
% We will use a patch object to fill the region between the outer and inner skin lines.
% The inner surface is defined by offsetting the outer surface by t/2 (simplified visualization).
y_inner_top = y_profile - thickness_profile / 2;
y_inner_bottom = -y_profile + thickness_profile / 2;

% Outer surface coordinates
x_outer = [x_chord, fliplr(x_chord)];
y_outer = [y_profile, fliplr(-y_profile)];

% 5. Create the plot using a series of 'Patch' objects
figure;
theme(gcf, "light");
hold on;
axis equal;

% Define the normalized color data (CData) based on thickness
% We normalize the thickness so that the color gradient spans the full range [0, 1]
t_min = min(thickness_profile);
t_max = max(thickness_profile);
normalized_t = (thickness_profile - t_min) / (t_max - t_min);

% Plot the skin section by section (patches) to enable the color gradient
for i = 1:(N - 1)
    % Define the vertices for a small trapezoidal section of the skin
    x_patch = [x_chord(i), x_chord(i+1), x_chord(i+1), x_chord(i)];
    
    % --- Top Skin Patch ---
    y_patch_top = [y_inner_top(i), y_inner_top(i+1), y_profile(i+1), y_profile(i)];
    % CData is the average normalized thickness for the color
    c_data_top = [normalized_t(i), normalized_t(i+1), normalized_t(i+1), normalized_t(i)];
    patch(x_patch, y_patch_top, c_data_top, 'EdgeColor', 'interp');

    % --- Bottom Skin Patch ---
    y_patch_bottom = [-y_profile(i), -y_profile(i+1), y_inner_bottom(i+1), y_inner_bottom(i)];
    c_data_bottom = [normalized_t(i), normalized_t(i+1), normalized_t(i+1), normalized_t(i)];
    patch(x_patch, y_patch_bottom, c_data_bottom, 'EdgeColor', 'interp');
end

% Set the colormap (e.g., jet or turbo for clear distinction)
colormap('jet');
c = colorbar;

% Labeling and Aesthetics
xlabel('Distance along Chord ($x$), m', 'Interpreter', 'latex');
ylabel('Height ($y$), m', 'Interpreter', 'latex');
title('2D NACA0012 Cross-Section with Optimal Skin Thickness Gradient');
c.Label.String = ['Skin Thickness, $t$ (m) | Start: ', num2str(t_max, '%.4e'), ' | End: ', num2str(t_min, '%.4e')];
c.Label.Interpreter = 'latex';

% Optional: Set limits for better viewing
y_max = max(y_profile) * 1.5;
ylim([-y_max, y_max]); 

hold off;

%% --- Objective Function Definition ---
function [f, maxStress, volume] = objectiveFunction(X, chord, dist_from_neutral_axis, N, targetStress, penaltyFactor, engFuncs)
    
    % Ensure thickness is non-negative
    if X(1) <= 0
        f = Inf;
        maxStress = Inf;
        volume = Inf;
        return;
    end
    
    startThickness = X(1);
    gradient = X(2);

    % 1. Evaluate the design (get stress and volume)
    [maxStress, volume] = evaluateDesign(startThickness, gradient, chord, dist_from_neutral_axis, N, engFuncs);

    % 2. Calculate the Stress Penalty
    stressRatio = (maxStress - targetStress) / targetStress;
    
    % Penalty is 0 if stress is met, otherwise it's a scaled squared error
    if stressRatio > 0
        penalty = penaltyFactor * stressRatio^2;
    else
        penalty = 0;
    end
    
    % 3. Calculate the objective value
    f = volume + penalty;
    
    fprintf('t0=%.4f, g=%.4f -> Volume: %.6e, MaxStress: %.3f MPa, Obj: %.6e\n', ...
            startThickness, gradient, volume, maxStress/1e6, f);
end


%% --- Helper Function to Evaluate Stress and Volume ---
function [maxStress, volume] = evaluateDesign(startThickness, gradient, chord, dist_from_neutral_axis, N, engFuncs)

    % 1. Compute I array
    I_array = engFuncs.secondMomentAreaArraySkinThicknessGradient(N, chord, startThickness, gradient);
    
    % Check for non-positive I (which can happen if thickness goes below zero)
    if any(I_array <= 0)
        maxStress = Inf;
        volume = Inf;
        return;
    end

    % 2. Compute Volume
    volume = engFuncs.findVolumeSkinThicknessGradient(chord, startThickness, gradient);

    % 3. Run beam bending solver
    % The beamBending function is assumed to be in your MATLAB path
    [~, M, ~, ~, ~] = beamBending(N, I_array);
    
    % 4. Compute Bending Stresses
    bendingStresses = engFuncs.bendingStress(M, dist_from_neutral_axis, I_array);

    % 5. Determine the maximum absolute stress
    maxStress = max(abs(bendingStresses));
end