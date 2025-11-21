% plot the profile
clc, clear all, close all;


engFuncs = makeEngFuncs;

N = 100;
chord = 0.15;
thickness_ratio = 0.12;
max_thickness = chord * thickness_ratio;
lin = linspace(0,chord, N);
yt = engFuncs.findAerofoilPositiveHeight(lin);

plot(lin, yt);
xlabel('Position along the chord (m)');
ylabel('Positive Height (m)');
title('Aerofoil Profile');
ylim([-0.075,0.075])
grid on;