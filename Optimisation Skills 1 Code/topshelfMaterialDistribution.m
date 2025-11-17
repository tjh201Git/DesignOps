% topshelf material distribution

clc, clear all, close all;

file = load('bestSecondMomentAreas.mat');
N = length(file.bestSecondMomentAreas);

engFuncs = makeEngFuncs;

shelfHeight = 0.001;

chord = 0.15;
heightRatio = 0.12;
maxPositiveHeight = chord*heightRatio/2;



% xPosMaxHeight = findXPosMaxHeight(engFuncs);



% now do binary search to find x pos start and end at certain height offset
% from neutral axis

% binary_search_tolerance = 0.00001;
% binary_search_max_iterations = 1000;

% Perform binary search to find the x positions corresponding to the shelf height
% xstart = findXPosAtHeightBinarySearch(engFuncs, shelfHeight, 0, xPosMaxHeight, binary_search_tolerance, binary_search_max_iterations)
% xend = findXPosAtHeightBinarySearch(engFuncs, shelfHeight, xPosMaxHeight, chord, binary_search_tolerance, binary_search_max_iterations)

% h1 = engFuncs.findAerofoilPositiveHeight(xpos_start)
% h2 = engFuncs.findAerofoilPositiveHeight(xpos_end)

% topShelfArea = findTopShelfArea(engFuncs, xstart, xend, shelfHeight)
% 
% topShelfSecondMomentArea = findTopShelfSecondMomentArea(xstart,xend, 1000, shelfHeight, engFuncs)

% maxsecondmomentarea = findTopShelfSecondMomentArea(0,chord, 0, 1000,  engFuncs)


% function secondMomentArea = findTopShelfSecondMomentArea(xstart, xend, shelfHeight, slices, engFuncs)
% 
%     dx = (xend - xstart) / slices;
%     dxoffset = dx/2;
%     xvals = xstart + dxoffset:dx:xend - dxoffset;
% 
%     aeroHeights = engFuncs.findAerofoilPositiveHeight(xvals);
%     rectHeights = aeroHeights - shelfHeight;
%     rectWidth = dx;
%     localSecondMomentAreas = engFuncs.secondMomentAreaRectangle(rectWidth, rectHeights);
% 
%     rectArea = rectWidth * rectHeights;
%     rectMidpointHeight = shelfHeight + rectHeights / 2
%     secondMomentAreaParts = engFuncs.applySecondMomentOffset(localSecondMomentAreas, rectArea, rectMidpointHeight)
% 
%     secondMomentArea = sum(secondMomentAreaParts) * 2 % mult by two for top and bottom of aerofoil
% 
% end

function I = findTopShelfSecondMomentArea(engFuncs, shelfHeight, xstart, xend, slices)

    dx = (xend - xstart)/slices;
    xvals = linspace(xstart + dx/2, xend - dx/2, slices);

    aeroHeights = engFuncs.findAerofoilPositiveHeight(xvals);
    rectHeights = aeroHeights - shelfHeight;

    % local rectangle area
    A = rectHeights .* dx;

    % local second moment about centroid
    I_local = engFuncs.secondMomentAreaRectangle(dx, rectHeights);

    % centroid position above neutral axis
    y_centroid = shelfHeight + rectHeights/2;

    % apply parallel-axis theorem
    I_total = I_local + A .* (y_centroid.^2);

    % double for symmetry
    I = 2 * sum(I_total);
end

function xPos = findXPosAtHeightBinarySearch(engFuncs, shelfHeight, xsearchstart, xsearchend, tolerance, iterations)

    xPos = binarySearch(@binarySearchHeightDiff, xsearchstart, xsearchend, tolerance, iterations);

    function height = binarySearchHeightDiff(x)
        aeroHeight = engFuncs.findAerofoilPositiveHeight(x)
        height = aeroHeight - shelfHeight;
    end
end

function area = findTopShelfArea(engFuncs, xstart, xend, shelfHeight)

    area = integral(@findHeightDifference, xstart, xend);

    function height = findHeightDifference(x)
        height = engFuncs.findAerofoilPositiveHeight(x) - shelfHeight;
    end
end

function xPosMaxHeight = findXPosMaxHeight(engFuncs)
    negHeight = @(x) -engFuncs.findAerofoilPositiveHeight(x);
    negGradient = @(x) -engFuncs.findAerofoilHeightDifferential(x);

    chord = 0.15;
    heightRatio = 0.12;
    maxPositiveHeight = chord * heightRatio;

    
    tolerance = 0.00001;
    max_iterations = 1000;
    alpha = 0.01;
    x0 = 0.3*chord;
    % gradien
    [xPosMaxHeight,~,~] = gradientDescent(negHeight, negGradient,x0,alpha, tolerance, max_iterations);
    height_validator = engFuncs.findAerofoilPositiveHeight(xPosMaxHeight);
    
    % if (abs(height_validator - maxPositiveHeight) < tolerance)
       % fprintf("Correct max height (%d) found at %d metres",maxPositiveHeight, xPosMaxHeight);
    % else
       % fprintf("Incorrect %d metres", xPosMaxHeight);
    % end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FIND xstart AND xend WHERE THE SHELF INTERSECTS THE AIRFOIL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xstart, xend] = findShelfSpan(engFuncs, shelfHeight, chord, xPosMaxHeight)

    % --- Max thickness location (approx or use your gradient descent) ----
    % xMax = 0.3 * chord;
    xMax = xPosMaxHeight;

    % Heights at the boundaries
    h_LE = engFuncs.findAerofoilPositiveHeight(0);
    h_max = engFuncs.findAerofoilPositiveHeight(xMax);
    h_TE = engFuncs.findAerofoilPositiveHeight(chord);

    % -------------------------------
    % LEFT INTERSECTION (xstart)
    % -------------------------------

    if shelfHeight <= h_LE + 1e-12
        % Shelf is at or below leading-edge height
        xstart = 0;

    elseif shelfHeight >= h_max - 1e-12
        % Shelf is above max thickness
        xstart = xMax;

    else
        % Shelf intersects between 0 and xMax
        xstart = binarySearch(@heightDiffLeft, 0, xMax, 1e-7, 200);
    end

    function y = heightDiffLeft(x)
        y = engFuncs.findAerofoilPositiveHeight(x) - shelfHeight;
    end



    % -------------------------------
    % RIGHT INTERSECTION (xend)
    % -------------------------------

    if shelfHeight <= h_TE + 1e-12
        % Shelf is below trailing edge height → intersection at TE
        xend = chord;

    elseif shelfHeight >= h_max - 1e-12
        % Shelf is above max thickness → only touches at xMax
        xend = xMax;

    else
        % Shelf intersects between xMax and chord
        xend = binarySearch(@heightDiffRight, xMax, chord, 1e-7, 200);
    end

    function y = heightDiffRight(x)
        y = engFuncs.findAerofoilPositiveHeight(x) - shelfHeight;
    end

end




function optimalShelfHeight = findShelfHeightForTargetSMA(engFuncs, idealSMA)

    % Geometry constants
    chord = 0.15;
    heightRatio = 0.12;
    maxPositiveHeight = chord * heightRatio;

    % Search bounds for shelf height
    h_min = 0;                % shelf fully at neutral axis
    h_max = maxPositiveHeight; % shelf at top surface

    xPosMaxHeight = findXPosMaxHeight(engFuncs);

    tol = 1e-6;
    maxIter = 200;

    % Outer binary search on shelf height
    optimalShelfHeight = binarySearch(@smaError, h_min, h_max, tol, maxIter);

    % The function we are trying to zero:
    function err = smaError(h)
        [xstart, xend] = findShelfSpan(engFuncs, h, chord, xPosMaxHeight);

        SMA = findTopShelfSecondMomentArea(engFuncs, h, xstart, xend, 800);
        err = SMA - idealSMA;   % want this = 0
    end
end

% optimalSMA = file.bestSecondMomentAreas(1)
% optimalShelfHeight = findShelfHeightForTargetSMA(engFuncs, optimalSMA)

% h0 = engFuncs.findAerofoilPositiveHeight(0)
% hc = engFuncs.findAerofoilPositiveHeight(0.15)

idealSMAs = file.bestSecondMomentAreas;
idealSMAsLength = length(idealSMAs);
idealShelfHeights = zeros(1,idealSMAsLength);

xPosMaxHeight = findXPosMaxHeight(engFuncs);

% lin = linspace(0,0.15,idealSMAsLength);

% plot(lin, idealShelfHeights);

fig = figure;
theme(fig, 'light');

hold on

lin = linspace(0,chord, idealSMAsLength);
yt = engFuncs.findAerofoilPositiveHeight(lin);



colors = parula(idealSMAsLength)

for i = 1:idealSMAsLength
    idealHeight = findShelfHeightForTargetSMA(engFuncs, idealSMAs(i));
    idealShelfHeights(i) = idealHeight;
    [xstart, xend] = findShelfSpan(engFuncs, idealHeight, chord, xPosMaxHeight);
    
plot([xstart, xend], [idealHeight, idealHeight], ...
         'Color', colors(i,:), ...
         'LineWidth', 1.0, ...
         'DisplayName', sprintf('Distance %d', lin(i)));

plot([xstart, xend], [-idealHeight, -idealHeight], ...
         'Color', colors(i,:), ...
         'LineWidth', 1.0, ...
         'DisplayName', sprintf('Distance %d', lin(i)));
end


plot(lin, yt, "Color", "black");
plot(lin, -yt, "Color", "black");

ylim([-0.075,0.075])

hold off

fig2 = figure;
theme(fig2, 'light');
hold on
grid on

colors = parula(idealSMAsLength);

% Scale SMA Index to 0 → 1.5 m
SMAspan = linspace(0, 1.5, idealSMAsLength);

%% -------------------------
% 1. Plot shelf lines in 3D
%% -------------------------

for i = 1:idealSMAsLength

    idealHeight = idealShelfHeights(i);           % use actual height
    [xstart, xend] = findShelfSpan(engFuncs, idealHeight, chord, xPosMaxHeight);

    % Positive shelf height
    plot3([xstart xend], [SMAspan(i) SMAspan(i)], [idealHeight idealHeight], ...
          'Color', colors(i,:), 'LineWidth', 1.5);

    % Negative shelf height
    plot3([xstart xend], [SMAspan(i) SMAspan(i)], [-idealHeight -idealHeight], ...
          'Color', colors(i,:), 'LineWidth', 1.5);
end

%% ----------------------------------------------------
% 2. Plot aerofoil shape (yt and -yt) across the SMA span
%% ----------------------------------------------------

[Xaf, Yaf] = meshgrid(lin, SMAspan);
Zupper = repmat(yt, idealSMAsLength, 1);
Zlower = -Zupper;

% Aerofoil upper surface
surf(Xaf, Yaf, Zupper, ...
     'FaceColor', [0.2 0.2 0.2], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.25);

% Aerofoil lower surface
surf(Xaf, Yaf, Zlower, ...
     'FaceColor', [0.2 0.2 0.2], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.25);

%% ----------------------------------
% Labels, limits, and view
%% ----------------------------------

xlabel('Chord Position (x)')
ylabel('Span (m)')
zlabel('Height (m)')
title('3D Shelf Lines + Aerofoil Shape Along Span 0 → 1.5 m')

colormap(parula(idealSMAsLength))
% zlim([-max(abs(idealShelfHeights)), max(abs(idealShelfHeights))]) % auto scale height
view(45, 25)
hold off


