function funcs = makeEngFuncs
    funcs.findVolumeSkinMethod = @findVolumeSkinMethod;
    funcs.findSecondMomentAreaSkinMethod = @findSecondMomentAreaSkinMethod;
    funcs.bendingStress = @bendingStress;
    funcs.secondMomentAreaArraySkinThicknessGradient = @secondMomentAreaArraySkinThicknessGradient;
    funcs.findAerofoilPositiveHeight = @findAerofoilPositiveHeight;
    funcs.findAerofoilHeightDifferential = @findAerofoilHeightDifferential;
    funcs.secondMomentAreaRectangle = @secondMomentAreaRectangle;
    funcs.applySecondMomentOffset = @applySecondMomentOffset;
    funcs.findAerofoilHeight = @findAerofoilHeight;
    funcs.findVolumeSkinThicknessGradient = @findVolumeSkinThicknessGradient;
    funcs.findAreaSkinThicknessGradient = @findAreaSkinThicknessGradient;

    funcs.approxCrossSectAreaSkinMethod = @approxCrossSectAreaSkinMethod;


    funcs.secondMomentAreaArraySkinThicknessQuadratic = @secondMomentAreaArraySkinThicknessQuadratic;
    funcs.findVolumeSkinThicknessQuadratic = @findVolumeSkinThicknessQuadratic;
    funcs.getThicknessAtDistanceFromQuadraticLine = @getThicknessAtDistanceFromQuadraticLine;
    funcs.findAreaAtDistanceSkinThicknessQuadratic = @findAreaAtDistanceSkinThicknessQuadratic;


    funcs.findSolidCrossSectionalArea = @findSolidCrossSectionalArea;
end
% Plot the aerofoil profile
% clc, clear all, close all
% 
% chord = 0.15;
% thickness = 1e-3;
% 
% n = 100;
% 
% x = linspace(0, chord, n);
% y = findAerofoilSurface(x, chord, thickness);
% yneg = -y;
% 
% fig = figure;
% theme(fig, "light");
% 
% hold on
% plot(x,y);
% plot(x,-y);
% y_axis = 5e-4;
% axis([0, chord, -y_axis, y_axis])
% 
% 
% hold off


% for NACA0012
function yt = findAerofoilPositiveHeight(x)

    % note that our aerofoil is symmetrical about the x axis

    chord = 0.15;
    thickness_ratio = 0.12;

    comp = 0.2969 * sqrt(x / chord);
    comp = comp - 0.1260 * (x / chord);
    comp = comp - 0.3516 * (x / chord).^2;
    comp = comp + 0.2843 * (x / chord).^3;
    comp = comp - 0.1015 * (x / chord).^4;

    yt = chord * 5 * thickness_ratio * comp;
   
end

function spurHeight = findAerofoilHeight(x)

    spurHeight = findAerofoilPositiveHeight(x) * 2;

end

function area = findSolidCrossSectionalArea(chord)
    area = integral(@findAerofoilHeight,0,chord);
end

function bending_stress = bendingStress(moment, dist_from_netural_axis, second_moment_area)
    bending_stress = moment .* dist_from_netural_axis ./ second_moment_area;
end

function differential = findAerofoilHeightDifferential(x)
    

chord = 0.15;
thickness_ratio = 0.12;

comp = 0.2969 ./ (sqrt(chord) * 2 * sqrt(x));
comp = comp - 0.1260 / chord;
comp = comp - 0.3516 * 2 * x / chord^2;
comp = comp + 0.2843 * 3 * x.^2 / chord^3;
comp = comp - 0.1015*4*x.^3 / chord^4;

differential = 5 * chord * thickness_ratio * comp;
end



% Approach 2 skin thickness (Skin method)
function area = approxCrossSectAreaSkinMethod(skin_thickness, chord)
    
    area = skin_thickness * integral(@integralPart, 0, chord) * 2;
    
    function out = integralPart(x)
        differential = findAerofoilHeightDifferential(x);
        out = sqrt(1+differential.^2);
    end

end

% To find volume we could have the same thickness down the length of the
% aerofoil, or we could have it in steps, or we could model it with a
% function. Since we want optimisation, i think we should do it in small
% chunks / steps
function V = findVolumeSkinMethod(chord, thickness, span)

    % num_thicknesses = length(thicknesses);
    % chunk_length = chord / num_thicknesses;
    % 
    % V = 0;
    % for i = 1:num_thicknesses
    %     V = V + approxCrossSectAreaSkinMethod(thicknesses(i), chord) * chunk_length;
    % end

    V = approxCrossSectAreaSkinMethod(thickness, chord) * span;

end


function secondMomentAreaArray = secondMomentAreaArraySimpleSkin(N, chord, thickness)
    secondMomentAreaArray = findSecondMomentAreaSkinMethod(chord, thickness) * ones(N);
end


% assumes single thickness per cross section
function I = findSecondMomentAreaSkinMethod(chord, thickness)
    
    numPanels = 100;
    panelLength = chord/numPanels;


    I = 0;
    for i = 0:(numPanels-1)
        x = (panelLength / 2) + (i * panelLength);
        height = findAerofoilPositiveHeight(x);
        differential = findAerofoilHeightDifferential(x);
        angle = atan(differential);
        I = I + secondMomentAreaPerSkinPanel(panelLength, thickness, angle, height);

    end
    % we have done it for the top side. Since it is symmetrical, just mult
    % by 2
    I = 2 * I; % Account for symmetry in the second moment of area
    
end


% L is length
% t is thickness
% angle is angle
% d is distance from centroid of rectangle to neutral axis
function Ix = secondMomentAreaPerSkinPanel(L, t, angle, d)
 Ix_ = L^3 * t * sin(angle)^2 / 12;
 A = t * L;

 Ix = Ix_ + A*d^2;
end


function thickness = getThicknessAtDistanceFromLinearLine(distance, startThickness, gradient)
    tempThickness = startThickness + gradient * distance;
    positiveThickness = max(0,tempThickness);
    thickness = min(positiveThickness, startThickness);
end

function area = findAreaSkinThicknessGradient(x, startThickness, gradient, chord)
        thickness = getThicknessAtDistanceFromLinearLine(x, startThickness, gradient);
        area = approxCrossSectAreaSkinMethod(thickness, chord);
end

function volume = findVolumeSkinThicknessGradient(chord, startThickness, gradient, span)
    
    volume = integral(@getAreaAtDistance, 0, span);

    function area = getAreaAtDistance(x)
        area = findAreaSkinThicknessGradient(x, startThickness, gradient, chord);
    end
end

function secondMomentAreaArray = secondMomentAreaArraySkinThicknessGradient(N, chord, startThickness, gradient)
    
    chunkDist = chord / N;
    % secondMomentAreaArray = zeros(N);

    for i = 1:N
        zpos = chunkDist * (i - 1);
        thickness = getThicknessAtDistanceFromLinearLine(zpos, startThickness, gradient);
        secondMomentAreaArray(i) = findSecondMomentAreaSkinMethod(chord, thickness);
    end
    
end

function secondMomentArea =  secondMomentAreaRectangle(width, height)
    secondMomentArea = width.*(height.^3) / 12;
end

function secondMomentArea = applySecondMomentOffset(localSecondMomentArea, area, distanceFromNeutralAxis)

    secondMomentArea = localSecondMomentArea + area .* distanceFromNeutralAxis.^2;
end








%% ---------------- Quadratic Thickness Helper Functions ----------------

function thickness = getThicknessAtDistanceFromQuadraticLine(distance, a, b, c)
    % Compute quadratic thickness at a given distance
    tempThickness = a + b * distance + c * distance.^2;
    % enforce non-negative thickness
    positiveThickness = max(0, tempThickness);
    % optional: limit thickness to root thickness (like in linear version)
    thickness = min(positiveThickness, a);
end

function area = findAreaAtDistanceSkinThicknessQuadratic(x, a, b, c, chord)
    thickness = getThicknessAtDistanceFromQuadraticLine(x, a, b, c);
    area = approxCrossSectAreaSkinMethod(thickness, chord);
end


function volume = findVolumeSkinThicknessQuadratic(chord, span, a, b, c)
    % Compute total skin volume by integrating cross-sectional area along chord
    volume = integral(@getAreaAtDistance, 0, span);

    function area = getAreaAtDistance(x)
        area = findAreaAtDistanceSkinThicknessQuadratic(x, a, b, c, chord);
    end
end


function secondMomentAreaArray = secondMomentAreaArraySkinThicknessQuadratic(N, chord, a, b, c)
    % Compute second moment of area array for N segments along chord
    chunkDist = chord / N;
    secondMomentAreaArray = zeros(1,N);  % preallocate

    for i = 1:N
        zpos = chunkDist * (i - 1);
        thickness = getThicknessAtDistanceFromQuadraticLine(zpos, a, b, c);
        secondMomentAreaArray(i) = findSecondMomentAreaSkinMethod(chord, thickness);
    end
end
