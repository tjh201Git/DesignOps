% Nelder mead optimised gradient skin thickness WITH BOX CONSTRAINTS

clc, clear all, close all

engFuncs = makeEngFuncs;

% --- Fixed Parameters ---
chord = 0.15;
span = 1.5;
dist_from_neutral_axis = 0.12 * chord;
N = 100;
targetStress = -30e6; % Pa
penaltyFactor = 1e7;

% --- VALID Initial Guess (must satisfy t0 + g >= 0) ---
X0 = [0.02, -0.005];   % SAFE FEASIBLE START

% --- BOX CONSTRAINTS ---
t0_min = 0.001;
t0_max = 0.050;

g_min = -0.05;   % more negative allowed
g_max = -0.000001; % always negative  

lb = [t0_min, g_min];
ub = [t0_max, g_max];

% --- Objective Function Handle ---
objective = @(X) objectiveFunction(X, chord, span, dist_from_neutral_axis, N, targetStress, penaltyFactor, engFuncs);

% --- Optimisation Options ---
options = optimset('Display','iter', 'TolX',1e-3, 'TolFun',1e-4);

% --- Run Bounded Nelder-Mead ---
[X_opt, fval] = fminsearchbnd(objective, X0, lb, ub, options);

% --- Extract Optimised Parameters ---
startThickness_opt = X_opt(1);
gradient_opt = X_opt(2);

% --- Recompute final design ---
[maxStress_opt, volume_opt] = evaluateDesign(startThickness_opt, gradient_opt, chord, span, dist_from_neutral_axis, N, engFuncs);

fprintf('\n*** Optimization Complete ***\n');
fprintf('Optimal Start Thickness: %.4f m\n', startThickness_opt);
fprintf('Optimal Gradient: %.4f\n', gradient_opt);
fprintf('Minimum Volume: %.6e m^3\n', volume_opt);
fprintf('Maximum Bending Stress: %.3f MPa\n', maxStress_opt / 1e6);
fprintf('Target Stress (magnitude): %.0f MPa\n', abs(targetStress) / 1e6);


%% --------------------- PLOTS ---------------------
figure;
theme(gcf, "light");

subplot(2,1,1);
dists = linspace(0,1,N);
spans = linspace(0,span,N);
thicknesses = startThickness_opt + spans * gradient_opt;

thicknesses = thicknesses * 1000;
plot(dists, thicknesses, 'LineWidth', 2);
xlabel('Normalised Radial position, r/R');
ylabel('Thickness, mm');
title("Optimised (Nelder-Mead) First Order Polynomial Thickness Over Aerofoil");


subplot(2,1,2)
I_array_opt = engFuncs.secondMomentAreaArraySkinThicknessGradient(N, chord, startThickness_opt, gradient_opt);
[~, M, ~, ~, ~] = beamBending(N, I_array_opt);
bendingStresses_opt = engFuncs.bendingStress(M, dist_from_neutral_axis, I_array_opt);

hold on;
plot(dists, abs(bendingStresses_opt)/1e6, 'b-', 'LineWidth', 2);
line([0,1], [abs(targetStress)/1e6, abs(targetStress)/1e6], 'Color','r','LineStyle','--');

title("Absolute Bending Stress Over Aerofoil Span")
xlabel('Normalised Radial position, r/R');
ylabel('Bending Stress, MPa');
legend('Computed |Stress|', 'Target Stress (30 MPa)', 'Location', 'best');
grid on;
hold off;



%% ---------------- Objective Function ----------------
function [f, maxStress, volume] = objectiveFunction(X, chord, span, dist_from_neutral_axis, N, targetStress, penaltyFactor, engFuncs)

    startThickness = X(1);
    gradient       = X(2);

    % enforce non-negative thickness everywhere
    t_min = min(startThickness, startThickness + gradient);
    if t_min <= 0
        f = Inf; maxStress = Inf; volume = Inf;
        return;
    end

    % Evaluate design
    [maxStress, volume] = evaluateDesign(startThickness, gradient, chord, span, dist_from_neutral_axis, N, engFuncs);

    % stress magnitude comparison
    targetMag = abs(targetStress);
    maxStressMag = abs(maxStress);

    if targetMag == 0
        stressRatio = 0;
    else
        stressRatio = (maxStressMag - targetMag) / targetMag;
    end

    if stressRatio > 0
        penalty = penaltyFactor * stressRatio^2;
    else
        penalty = 0;
    end

    f = volume + penalty;

    fprintf('t0=%.4f, g=%.4f -> Volume: %.6e, MaxStress: %.3f MPa, Obj: %.6e\n', ...
        startThickness, gradient, volume, maxStressMag/1e6, f);
end


%% ---------------- Evaluate Design ----------------
function [maxStress, volume] = evaluateDesign(startThickness, gradient, chord, span, dist_from_neutral_axis, N, engFuncs)

    % enforce t >= 0 everywhere
    if startThickness + gradient*span < 0
        maxStress = Inf;
        volume = Inf;
        return
    end

    I_array = engFuncs.secondMomentAreaArraySkinThicknessGradient(N, chord, startThickness, gradient);

    if any(I_array <= 0)
        maxStress = Inf;
        volume = Inf;
        return;
    end

    volume = engFuncs.findVolumeSkinThicknessGradient(chord, startThickness, gradient, span);

    [~, M, ~, ~, ~] = beamBending(N, I_array);

    bendingStresses = engFuncs.bendingStress(M, dist_from_neutral_axis, I_array);

    maxStress = max(abs(bendingStresses));
end
