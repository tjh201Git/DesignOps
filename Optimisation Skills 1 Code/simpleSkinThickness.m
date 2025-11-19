% Simple same thickness the entire way through
clc, clear all, close all

engFuncs = makeEngFuncs;

youngs_modulus = 2e9;
yield_strength = 30e6;

% N = 100;
% chord = 0.15;
% thickness = 0.03;



% --- Parameters
chord = 0.15;
dist_from_neutral_axis = 0.12*chord;

N = 100;
targetStress = 30e6; % Pa (30 MPa)

% Search bounds for thickness (m)
t_low = 0.00005;  
t_high = 0.01;   

% Run binary search
[t_opt, history] = findOptimalThicknessBinary(t_low, t_high, chord, engFuncs, N, targetStress);

% fprintf('\nFinal optimal thickness = %.6e m (%.3f mm)\n', t_opt, t_opt*1e3);
% fprintf('\nFinal ')

iters_num = length(history.t);
iters = 1:iters_num;

fig = figure;
theme(fig, "light");

% --- Left Y-Axis (Thickness) ---
yyaxis left;
plot(iters, history.t, 'o-'); % Plot Thickness
xlabel('Iterations, N');
ylabel('Thickness, t');

% --- Right Y-Axis (Volume) ---
yyaxis right;
plot(iters, history.volume, '*-'); % Plot Volume
% Set the right y-axis label to 'Volume, m^3' using the LaTeX interpreter
ylabel('Volume, m$^3$', 'Interpreter', 'latex');
title("Binary Search For Smallest Uniform Thickness")
% Determine the min/max limits for the volume data
% V_min = min(history.volume);
% V_max = max(history.volume);

% Apply these limits to the right y-axis scale.
% This automatically maps the left y-axis's plot range to the
% corresponding volume range on the right.
% ylim([V_min, V_max]);


finalThickness = history.t(end);
% Compute volume using your engineering function
finalVolume = engFuncs.findVolumeSkinMethod(chord, finalThickness);

% Compute second moment of area
final_I = engFuncs.findSecondMomentAreaSkinMethod(chord, finalThickness);
I_array = final_I .* ones(1, N);

% Run beam bending solver
[~, M, ~, ~, ~] = beamBending(N, I_array);

% Distance from neutral axis 
y = 0.12 * chord;

% Compute bending stresses
finalStresses = engFuncs.bendingStress(M, y, I_array);
finalStresses = finalStresses / 1e6;

lin = linspace(0,1,N);

myColourMap = flipud(jet(N));


fig2 = figure;
theme(fig2, "light");
% Plot points and color them based on their index (which is linear)
hold on
plot(lin, finalStresses, 'black');
scatter(lin, finalStresses, 30, lin, 'filled'); % 30 is the marker size, lin maps the color
% colorbar; % Display the color bar

% Apply your colormap and labels
colormap(flipud(jet(N)));
title("Bending Stress Along Aerofoil, Uniform Thickness 2.382mm")
xlabel('Normalised Radial position, r/R');
ylabel('Bending Stress, MPa');
% subplot(3,1,3);
% plot(history.volume, history.t, 'ko-', 'Color', 'b');
% xlabel('volume');
% ylabel('thickness');

hold off

function [t_opt, history] = findOptimalThicknessBinary(t_low, t_high, chord, engFuncs, beamN, targetStress, tol)
% findOptimalThicknessBinary
% Uses binary search to find the minimum skin thickness such that
% max bending stress <= targetStress.
%
% INPUTS:
%   t_low        - lower bound (m)
%   t_high       - upper bound (m)
%   chord        - chord length (m)
%   engFuncs     - struct returned by makeEngFuncs
%   beamN        - number of nodes for beamBending
%   targetStress - allowable maximum bending stress (Pa)
%   tol          - tolerance on thickness convergence (m)
%
% OUTPUTS:
%   t_opt   - optimal thickness (m)
%   history - struct containing all iteration data

    if nargin < 7
        tol = 1e-5; % default thickness tolerance (m)
    end

    maxIter = 100;
    iter = 0;

    history.t = [];
    history.maxStress = [];
    history.volume = [];

    while (t_high - t_low) > tol && iter < maxIter
        iter = iter + 1;
        t_mid = 0.5 * (t_low + t_high);

        [maxStress_mid, volume_mid] = evaluateThickness(t_mid, chord, engFuncs, beamN);

        history.t(end+1) = t_mid;
        history.maxStress(end+1) = maxStress_mid;
        history.volume(end+1) = volume_mid;

        if maxStress_mid > targetStress
            % too thin -> stresses too high -> increase thickness
            t_low = t_mid;
        else
            % safe (stress below target) -> try thinner
            t_high = t_mid;
        end

        fprintf("Iter %2d: t = %.6e m (%.3f mm), maxStress = %.3e Pa, volume = %.3e M^3 \n",  ...
                iter, t_mid, t_mid*1e3, maxStress_mid, volume_mid);
    end

    % final result
    t_opt = t_high;
    volume_opt = history.volume(end);
    fprintf("Optimal thickness found: %.6e m (%.3f mm), Volume = %.3e M^3 \n", t_opt, t_opt*1e3, volume_opt);

end


%% Helper function to compute stress and volume for given thickness
function [maxStress, volume] = evaluateThickness(t, chord, engFuncs, N)

    % Compute volume using your engineering function
    volume = engFuncs.findVolumeSkinMethod(chord, t);

    % Compute second moment of area
    I = engFuncs.findSecondMomentAreaSkinMethod(chord, t)
    I_array = I * ones(1, N);

    % Run beam bending solver
    [~, M, ~, ~, ~] = beamBending(N, I_array);

    % Distance from neutral axis 
    y = 0.12 * chord;

    % Compute bending stresses
    stresses = engFuncs.bendingStress(M, y, I_array);

    % Return the maximum absolute bending stress
    maxStress = max(abs(stresses));
end
