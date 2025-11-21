% calculate solid volume

engFuncs = makeEngFuncs;

chord = 0.15;
span = 1.5;

area = engFuncs.findSolidCrossSectionalArea(chord);
volume = area*span;

fprintf("Area, %.10f\n", area);
fprintf("Volume: %.10f\n", volume);
    