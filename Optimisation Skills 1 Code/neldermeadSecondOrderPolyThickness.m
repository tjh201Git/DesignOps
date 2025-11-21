% ---------------------------------------------------------
% Second-Order (Quadratic) Thickness Optimisation
% **Updated to Operate Over PHYSICAL Span x in [0, span]**
% ---------------------------------------------------------
clc; clear all; close all;
% NOTE: Assuming 'makeEngFuncs' and 'beamBendingFast' are available.
engFuncs = makeEngFuncs; 
% --- Fixed Parameters ---
chord = 0.15;
span = 1.5;
dist_from_neutral_axis = 0.12 * chord;
N = 100;
targetStress = -30e6; % Pa
penaltyFactor = 1e7;
% --- Target Tip Thickness ---
t_tip_target = 2e-6; % 0.02 mm (Change this to whatever minimum you want)
% --- VALID Initial Guess (Only a and b) ---
% Start with a reasonable root thickness (a) and slope (b)
X0 = [0.02, -0.015]; 
% --- BOX CONSTRAINTS ---
a_min = 0.001;  a_max = 0.05;
b_min = -0.10;  b_max = 0.00; % b must be negative to taper down
lb = [a_min, b_min];
ub = [a_max, b_max];
% Update the handle to pass t_tip_target
objective = @(X) objectiveFunctionQuad(X, t_tip_target, chord, span, dist_from_neutral_axis, N, targetStress, penaltyFactor, engFuncs);
% --- Optimisation Options ---
options = optimset('Display','iter', 'TolX',1e-3, 'TolFun',1e-4, 'MaxIter', 300);
% Run Optimization
[X_opt, fval] = fminsearchbnd(objective, X0, lb, ub, options);
% --- Extract and Calculate c ---
a_opt = X_opt(1);
b_opt = X_opt(2);
% C is calculated for PHYSICAL span x=[0, span]: t(span) = a + b*span + c*span^2 = t_tip_target
c_opt = (t_tip_target - a_opt - b_opt * span) / span^2; % <-- FIXED for PHYSICAL span
% --- Recompute final design ---
[maxStress_opt, volume_opt] = evaluateDesignQuad(a_opt, b_opt, c_opt, chord, span, dist_from_neutral_axis, N, engFuncs);
fprintf('\n*** Optimization Complete ***\n');
fprintf('Optimal a (root thickness): %.10f m\n', a_opt);
fprintf('Optimal b (linear term, PHYSICAL): %.10f\n', b_opt);
fprintf('Optimal c (quadratic term, PHYSICAL): %.10f\n', c_opt);
fprintf('Minimum Volume: %.10f m^3\n', volume_opt);
fprintf('Maximum Bending Stress: %.10f MPa\n', maxStress_opt / 1e6);
fprintf('Target Stress (magnitude): %.0f MPa\n', abs(targetStress) / 1e6);
%% --------------------- PLOTS ---------------------
fig = figure;
theme(fig, "light"); 
subplot(2,1,1);
spanDists = linspace(0,span,N); % Physical distance
thicknesses = a_opt + b_opt*spanDists + c_opt*spanDists.^2; % Uses physical distance coefficients
thicknesses = thicknesses * 1000;
linDists = spanDists / span; % Normalized distance for plotting
plot(linDists, thicknesses, 'LineWidth', 2);
xlabel('Normalised Radial Position, r/R');
ylabel('Thickness, mm');
title("Optimised (Nelder-Mead) Second Order Polynomial Thickness Over Aerofoil");
subplot(2,1,2);
% Call to engFuncs uses only the coefficients (a, b, c)
I_array_opt = engFuncs.secondMomentAreaArraySkinThicknessQuadratic(N, chord, a_opt, b_opt, c_opt); 
[~, M, ~, ~, ~] = beamBending(N, I_array_opt);
bendingStresses_opt = engFuncs.bendingStress(M, dist_from_neutral_axis, I_array_opt);
hold on;
plot(linDists, abs(bendingStresses_opt)/1e6, 'b-', 'LineWidth', 2);
line([0,1], [abs(targetStress)/1e6, abs(targetStress)/1e6], 'Color','r','LineStyle','--');
title("Absolute Bending Stress Over Aerofoil Span");
xlabel('Normalised Radial Position, r/R');
ylabel('Bending Stress, MPa');
legend('Computed |Stress|', 'Target Stress', 'Location', 'best');
grid on;
hold off;
%% ---------------- Objective Function ----------------
function [f, maxStress, volume] = objectiveFunctionQuad(X, t_tip_target, chord, span, dist_from_neutral_axis, N, targetStress, penaltyFactor, engFuncs)
    a = X(1); 
    b = X(2);
    
    % --- FORCE TIP THICKNESS ---
    % Calculate c so that at x=span, t = t_tip_target
    c = (t_tip_target - a - b * span) / span^2; % <-- FIXED for PHYSICAL span
    
    % --- 1. Monotonicity Check (over physical span x=[0,span]) ---
    % The derivative is t'(x) = 2*c*x + b. 
    deriv_0 = b;              % Derivative at x=0
    deriv_span = b + 2 * c * span;  % Derivative at x=span
    
    if deriv_0 > 0 || deriv_span > 0
        % If derivative is positive, thickness is increasing somewhere.
        f = Inf; maxStress = Inf; volume = Inf;
        return;
    end
    % --- 2. Convexity/Positivity Check (over physical span x=[0,span]) ---
    x_check = 0.5 * span; % Check at physical midpoint
    if (a + b*x_check + c*x_check^2) <= 0
        f = Inf; maxStress = Inf; volume = Inf;
        return;
    end
    % --- Evaluate design ---
    [maxStress, volume] = evaluateDesignQuad(a, b, c, chord, span, dist_from_neutral_axis, N, engFuncs);
    
    % --- Stress penalty ---
    targetMag = abs(targetStress);
    maxStressMag = abs(maxStress);
    
    if maxStressMag > targetMag
        stressRatio = (maxStressMag - targetMag)/targetMag;
        penalty = penaltyFactor * stressRatio^2;
    else
        penalty = 0;
    end
    
    f = volume + penalty;
end
%% ---------------- Evaluate Design ----------------
function [maxStress, volume] = evaluateDesignQuad(a, b, c, chord, span, dist_from_neutral_axis, N, engFuncs)
    
    % --- Evaluate thickness at critical points (over physical span x=[0,span]) ---
    x = linspace(0,span,N); 
    
    % The thickness equation uses physical coefficients (a, b, c) on physical distance (x)
    t_array = a + b*x + c*x.^2;
    
    % Ensure thickness never negative
    if any(t_array <= 0)
        maxStress = Inf; volume = Inf;
        return;
    end
    
    % Ensure monotonic decreasing (Derivative check over physical span)
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
    volume = engFuncs.findVolumeSkinThicknessQuadratic(chord, span, a, b, c);
    % --- Beam bending and stresses ---
    [~, M, ~, ~, ~] = beamBending(N, I_array);
    bendingStresses = engFuncs.bendingStress(M, dist_from_neutral_axis, I_array);
    maxStress = max(abs(bendingStresses));
end