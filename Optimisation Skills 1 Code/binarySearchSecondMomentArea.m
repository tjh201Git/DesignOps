% find exact second moment area at points along wing
clc, clear all, close all

engFuncs = makeEngFuncs;

youngs_modulus = 2e9;
yield_strength = 30e6; %MPa
chord = 0.15;
radius = 1.5;
dist_from_neutral_axis = 0.12*chord;

I_low = 1e-10;
I_high = 1e-7;
I_tol = 1e-12;

N = 100;

history = findOptimalThicknessBinary(I_low, I_high, chord, N, yield_strength, I_tol, engFuncs);
history.stresses = history.stresses ./ 1e6;

bestSecondMomentAreas = history.secondMomentAreas(end,:);
save('bestSecondMomentAreas.mat','bestSecondMomentAreas');


fig = figure;
theme(fig, "light");

% x-axis = spanwise or chordwise locations
lin = linspace(0, radius, N);

subplot(2,1,1);
hold on

% number of iterations stored
numIters = size(history.stresses, 1);

% plot each iteration with a color gradient
colors = parula(numIters);

for k = 1:numIters
    plot(lin, history.stresses(k,:), ...
         'Color', colors(k,:), ...
         'LineWidth', 1.0, ...
         'DisplayName', sprintf('Iter %d', k));
end

% highlight the last iteration (final result)
plot(lin, history.stresses(end,:), ...
     'k', 'LineWidth', 2.0, ...
     'DisplayName', 'Final Converged');

xlabel("Spanwise Position (m)");
ylabel("Bending Stress (MPa)");
title("Stress Distribution Convergence using Binary Search");

legend("show", "Location", "bestoutside");
grid on

hold off

subplot(2,1,2);
hold on

% number of iterations stored
numIters = size(history.stresses, 1);

% plot each iteration with a color gradient
colors = parula(numIters);

for k = 1:numIters
    plot(lin, history.secondMomentAreas(k,:), ...
         'Color', colors(k,:), ...
         'LineWidth', 1.0, ...
         'DisplayName', sprintf('Iter %d', k));
end

% highlight the last iteration (final result)
plot(lin, history.secondMomentAreas(end,:), ...
     'k', 'LineWidth', 2.0, ...
     'DisplayName', 'Final Converged');

xlabel("Spanwise Position (m)");
ylabel('Second Moment of Area (m^{4})')
title("Second Moment of Area Convergence using Binary Search");

legend("show", "Location", "bestoutside");
grid on

hold off


function history = findOptimalThicknessBinary(I_low, I_high, chord, N, targetStress, tol, engFuncs)
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

    maxIter = 100;
    iter = 0;

    % history.thicknesses = [];
    history.stresses = [];
    history.secondMomentAreas = [];
    % history.volume = [];

    I_highs = I_high * ones(1, N);
    I_lows = I_low * ones(1,N);

    % calculate max diff between high and low guesses
    I_diffs = I_highs - I_lows;
    I_max_diff = max(I_diffs);

    while I_max_diff > tol && iter < maxIter

        fprintf("Iteration %i\n", iter);



        iter = iter + 1;
        I_mids = 0.5 * (I_lows + I_highs);

        stresses = evaluateThicknesses(I_mids, chord);

        % % history.t(end+1) = t_mid;
        % history.stresses(end+1) = stresses;
        % % history.volume(end+1) = volume_mid;
        % history.secondMomentAreas(end+1) = I_mids;

        history.stresses(iter, :) = stresses(:).';          % store row
        history.secondMomentAreas(iter, :) = I_mids(:).';   % store row


        for i = 1:length(I_mids)
            if abs(stresses(i)) > targetStress
                % too thin -> stresses too high -> increase thickness
                I_lows(i) = I_mids(i);
            else
                % safe (stress below target) -> try thinner
                I_highs(i) = I_mids(i);
            end
        end
        % fprintf("Iter %2d: t = %.6e m (%.3f mm), maxStress = %.3e Pa, volume = %.3e M^3 \n",  ...
        %         iter, t_mid, t_mid*1e3, maxStress_mid, volume_mid);


        I_diffs = I_highs - I_lows;
        I_max_diff = max(I_diffs);
    end

    % final result
    % t_opt = t_high;
    % volume_opt = history.volume(end);
    % fprintf("Optimal thickness found: %.6e m (%.3f mm), Volume = %.3e M^3 \n", t_opt, t_opt*1e3, volume_opt);

end




%% Helper function to compute stress and volume for given thickness
function stresses = evaluateThicknesses(I_array, chord)

    N = length(I_array);

    % Compute volume using your engineering function
    % volume = engFuncs.findVolumeSkinMethod(chord, t);

    % Compute second moment of area
    % I = engFuncs.findSecondMomentAreaSkinMethod(chord, t);
    % I_array = I * ones(N, 1);

    % Run beam bending solvers
    [~, M, ~, ~, ~] = beamBending(N, I_array);

    % Distance from neutral axis 
    y = 0.12 * chord;

    % Compute bending stresses
    % stresses = engFuncs.bendingStress(M, y, I_array);
    stresses = M .* y ./ I_array;

    % Return the maximum absolute bending stress
    % maxStress = max(abs(stresses));
end
