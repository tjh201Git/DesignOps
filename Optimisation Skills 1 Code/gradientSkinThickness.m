% gradient skin thickness

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

startThickness = 0.05;
% gradient = -0.3;
gradient = -0.25;

I_array = engFuncs.secondMomentAreaArraySkinThicknessGradient(N, chord, startThickness, gradient);

dists = linspace(0, chord, N);

%submit to beamBending:
% number of nodes for numerical solve, N
% vector of second moment of area values across the span, I 
% (I must have N elements)
[delta, M, V, L, P]  = beamBending(N,I_array); 
%beamBending function returns:
%beam deflection, delta
%bending moment distribution, M
%shear force distriubition, V
%blade length is returned, L
%distributed aerodynamic loading, P

bendingStresses = engFuncs.bendingStress(M, dist_from_neutral_axis, I_array);





fig = figure;
theme(fig, "light");

% plot(dists, I, 'ko-');
plot(dists, bendingStresses);
xlabel('distance, m'); %set labels
ylabel('Bending Stress, MPa');




% subplot(3,1,1);
% plot(iters, history.t,'ko-', 'Color', 'g'); %plot the deflection (delta) at each node location across the blade
% xlabel('iterations'); %set labels
% ylabel('thickness');
% 
% subplot(3,1,2);
% plot(history.maxStress, history.t, 'ko-', 'Color', 'r');
% xlabel('Bending Stress, MPa');
% ylabel('thickness');
% 
% subplot(3,1,3);
% plot(history.volume, history.t, 'ko-', 'Color', 'b');
% xlabel('volume');
% ylabel('thickness');
% 
% function [t_opt, history] = findOptimalThicknessBinary(t_low, t_high, chord, engFuncs, beamN, targetStress, tol)
% % findOptimalThicknessBinary
% % Uses binary search to find the minimum skin thickness such that
% % max bending stress <= targetStress.
% %
% % INPUTS:
% %   t_low        - lower bound (m)
% %   t_high       - upper bound (m)
% %   chord        - chord length (m)
% %   engFuncs     - struct returned by makeEngFuncs
% %   beamN        - number of nodes for beamBending
% %   targetStress - allowable maximum bending stress (Pa)
% %   tol          - tolerance on thickness convergence (m)
% %
% % OUTPUTS:
% %   t_opt   - optimal thickness (m)
% %   history - struct containing all iteration data
% 
%     if nargin < 7
%         tol = 1e-5; % default thickness tolerance (m)
%     end
% 
%     maxIter = 100;
%     iter = 0;
% 
%     history.t = [];
%     history.maxStress = [];
%     history.volume = [];
% 
%     while (t_high - t_low) > tol && iter < maxIter
%         iter = iter + 1;
%         t_mid = 0.5 * (t_low + t_high);
% 
%         [maxStress_mid, volume_mid] = evaluateThickness(t_mid, chord, engFuncs, beamN);
% 
%         history.t(end+1) = t_mid;
%         history.maxStress(end+1) = maxStress_mid;
%         history.volume(end+1) = volume_mid;
% 
%         if maxStress_mid > targetStress
%             % too thin -> stresses too high -> increase thickness
%             t_low = t_mid;
%         else
%             % safe (stress below target) -> try thinner
%             t_high = t_mid;
%         end
% 
%         fprintf("Iter %2d: t = %.6e m (%.3f mm), maxStress = %.3e Pa, volume = %.3e M^3 \n",  ...
%                 iter, t_mid, t_mid*1e3, maxStress_mid, volume_mid);
%     end
% 
%     % final result
%     t_opt = t_high;
%     volume_opt = history.volume(end);
%     fprintf("Optimal thickness found: %.6e m (%.3f mm), Volume = %.3e M^3 \n", t_opt, t_opt*1e3, volume_opt);
% 
% end
% 
% 
% %% Helper function to compute stress and volume for given thickness
% function [maxStress, volume] = evaluateThickness(t, chord, engFuncs, N)
% 
%     % Compute volume using your engineering function
%     volume = engFuncs.findVolumeSkinMethod(chord, t);
% 
%     % Compute second moment of area
%     I = engFuncs.findSecondMomentAreaSkinMethod(chord, t);
%     I_array = I * ones(N, 1);
% 
%     % Run beam bending solver
%     [~, M, ~, ~, ~] = beamBending(N, I_array);
% 
%     % Distance from neutral axis 
%     y = 0.12 * chord;
% 
%     % Compute bending stresses
%     stresses = engFuncs.bendingStress(M, y, I);
% 
%     % Return the maximum absolute bending stress
%     maxStress = max(abs(stresses));
% end
