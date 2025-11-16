% max second moment area

engFuncs = makeEngFuncs;

chord = 0.15;
slices = 10000;
width = chord / slices;
x = linspace(0,chord, slices);
heights = engFuncs.findAerofoilHeight(x);
secondMomentAreas = engFuncs.secondMomentAreaRectangle(width, heights);
totalSMA = sum(secondMomentAreas)