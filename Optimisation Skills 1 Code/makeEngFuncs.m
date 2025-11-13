function funs = makeEngFuncs
    funs.findVolumeSkinMethodChunks = @findVolumeSkinMethodChunks;
    funs.findSecondMomentAreaSkinMethod = @findSecondMomentAreaSkinMethod;
    funs.bendingStress = @bendingStress;
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
function yt = findAerofoilPositiveHeight(x, chord, thickness)

    % note that our aerofoil is symmetrical about the x axis

    comp = 0.2969 * sqrt(x / chord);
    comp = comp - 0.1260 * (x / chord);
    comp = comp - 0.3516 * (x / chord).^2;
    comp = comp + 0.2843 * (x / chord).^3;
    comp = comp - 0.1015 * (x / chord).^4;

    yt = chord * 5 * thickness * comp;
   
end

function spurHeight = findAerofoilHeight(x, chord, thickness)

    spurHeight = findAerofoilPositiveHeight(x, chord, thickness) * 2;

end

function area = findCrossSectionalArea(chord, thickness)


    area = integral(@parameterisedFindAeroHeight,0,chord);
    
    % Use a nested function so we can use the matlab integral function with
    % our custom variables (chord and thickness)
    function height = parameterisedFindAeroHeight(x)
        height = 2 * findAerofoilHeight(x, chord, thickness);
    end

end

function bending_stress = bendingStress(moment, dist_from_netural_axis, second_moment_area)
    bending_stress = moment * dist_from_netural_axis ./ second_moment_area;
end

function differential = findAerofoilHeightDifferential(x, chord, thickness)
    
comp = 0.2969 ./ (sqrt(chord) * 2 * sqrt(x));
comp = comp - 0.1260 / chord;
comp = comp - 0.3516 * 2 * x / chord^2;
comp = comp + 0.2843 * 3 * x.^2 / chord^3;
comp = comp - 0.1015*4*x.^3 / chord^4;

differential = 5 * chord * thickness * comp;
end



% area = findCrossSectionalArea(0.15, 1e-3)



% Approach 2 skin thickness (Skin method)
function area = approxCrossSectAreaSkinMethod(skin_thickness, chord)
    
    area = skin_thickness * integral(@integralPart, 0, chord) * 2;
    
    function out = integralPart(x)
        differential = findAerofoilHeightDifferential(x, chord, skin_thickness);
        out = sqrt(1+differential.^2);
    end

end

% To find volume we could have the same thickness down the length of the
% aerofoil, or we could have it in steps, or we could model it with a
% function. Since we want optimisation, i think we should do it in small
% chunks / steps
function V = findVolumeSkinMethodChunks(chord, thicknesses)

    num_thicknesses = length(thicknesses);
    chunk_length = chord / num_thicknesses;

    V = 0;
    for i = 1:num_thicknesses
        V = V + approxCrossSectAreaSkinMethod(thicknesses(i), chord) * chunk_length;
    end

end


% assumes single thickness per cross section
function I = findSecondMomentAreaSkinMethod(chord, thickness)
    
    numPanels = 100;
    panelLength = chord/numPanels;


    I = 0;
    for i = 0:(numPanels-1)
        x = (panelLength / 2) + (i * panelLength);
        height = findAerofoilPositiveHeight(x, chord, thickness);
        differential = findAerofoilHeightDifferential(x, chord, thickness);
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


