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
gradient = -0.3;
% gradient = -0.25;

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

fminsearch()

