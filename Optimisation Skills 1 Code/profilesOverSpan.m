% plot different thickness profiles over aerofoil span

engFuncs = makeEngFuncs;

% mm
uniformThickness = 2.382;

firstOrderT0 = 2.3812;
firstOrderGradient = -1.59;

secondOrderA = 2.3882856;
secondOrderB = -3.18013;
secondOrderC = 1.0595153;

% *** Optimization Complete ***
% Optimal a (root thickness): 0.0023882856 m
% Optimal b (linear term, PHYSICAL): -0.0031801299
% Optimal c (quadratic term, PHYSICAL): 0.0010595153
% Minimum Volume: 0.0003660707 m^3
% Maximum Bending Stress: 29.9107359608 MPa
% Target Stress (magnitude): 30 MPa

chord = 0.15;
span = 1.5;

N = 100;

% Define the span of the aerofoil
spanArray = linspace(0, span, N);
linArray = linspace(0, 1, N);

% Calculate thickness profiles
uniformThicknessProfile = uniformThickness * ones(size(spanArray));
firstOrderThickness = firstOrderT0 + firstOrderGradient * spanArray;
secondOrderThickness = secondOrderA  + secondOrderB * spanArray + secondOrderC * spanArray.^2;

fig = figure;
theme(fig, 'light');

hold on

% Plot the thickness profiles
plot(linArray, uniformThicknessProfile, '-', 'DisplayName', 'Uniform Thickness');
plot(linArray, firstOrderThickness, '-', 'DisplayName', 'First Order Thickness');
plot(linArray, secondOrderThickness, '-', 'DisplayName', 'Second Order Thickness');

ylim([0, 2.5]);

title("Comparison of Different Structural Skin Thickness Profiles Over Aerofoil")
% colormap(winter);

ylabel("Structural Skin Thickness, mm")
xlabel("Normalised Radial Position, r/R");
legend

hold off


% fig2 = figure;
% theme(fig2, "light");
% 
% hold on
% 
% uniformArea = engFuncs.approxCrossSectAreaSkinMethod(uniformThickness, chord);
% firstOrderArea = engFuncs.findAreaSkinThicknessGradient(linArray, firstOrderT0, firstOrderGradient, chord);
% secondOrderArea = engFuncs.findAreaAtDistanceSkinThicknessQuadratic(linArray,secondOrderA,secondOrderB,secondOrderC,chord);
% 
% plot([0,1], [uniformArea,uniformArea], '--', 'DisplayName', 'Uniform Area');
% plot(linArray, firstOrderArea, '--', 'DisplayName', 'First Order Area');
% plot(linArray, secondOrderArea, '--', 'DisplayName', 'Second Order Area');
% legend;
% 
% 
% xlabel("Normalised Radial Position, r/R");
% ylabel("Area, m^2")
% 
% hold off