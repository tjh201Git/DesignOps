clc, clear all, close all

youngs_modulus = 2e9;
yield_strength = 30e6;

N = 100;
chord = 0.15;
dist_from_neutral_axis = 0.12*chord;
thicknesses = [0.0325];


% aim is to maximise second moment of area, whilst minimising volume

engFuncs = makeEngFuncs;

volume = engFuncs.findVolumeSkinMethodChunks(chord, thicknesses)

I = engFuncs.findSecondMomentAreaSkinMethod(chord, thicknesses(1));
I_array = I * ones(N);


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

%define the spanwise position of every node (in meters)
nodeLocations = linspace(0,L,N);

% Create a figure
fig = figure;

% Set the theme to dark
theme(fig, "light");

%plot the results
subplot(2,1,1) %tells MATLAB you want 2 subplots in one figure (check the docs for more info) 
plot(nodeLocations,delta,'ko-') %plot the deflection (delta) at each node location across the blade
xlabel('distance, m') %set labels
ylabel('deflection, mm')

subplot(2,1,2) %plot on the second subplot
bendingStress = (M*dist_from_neutral_axis./I)/(1e6);
plot(nodeLocations, bendingStress,'ko-') %plot the bending stress (M*y / I) at every node location
xlabel('distance, m') %set labels
ylabel('Bending Stress, MPa')