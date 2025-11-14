%MATLAB BLACKBOX OPTIMIZER
%=============================================================
%Steps
%
% 1) first samples the function using the Star discrepency model
% 2) Plots the result
% 3) Chose a cluster to begin Genetic Algorithm sampling
% 5) repeat until the highest MCSI is found
%====================================================================

minAdvert = -5;
maxAdvert = 5;
minTarget = -5;
maxTarget = 5;
minSkew = -5;
maxSkew = 5;

%===========Star Discrepency===================

Nodes  = 20;
Dimensions  = 3; 
               
% We want to generate many candidate point sets in normalised 3D space
% and measure the star discrepancy of each one. The goal is to find the
% set of points that is the most evenly space, which is the lowest 
% discrepency, to use as our initial sampling design.

numSets = 200; %Try 200 different designs
bestSD  = inf; %Start with infinity then replace with smaller and smaller discrepency values
bestXn  = []; %placeholder for best design

for k = 1:numSets %For each set

    Xn_cand = rand(Nodes, Dimensions); %Generate a random N×D point set in the hypercube
    sD_cand = starD(Xn_cand); %Find the star disrepency

    if sD_cand < bestSD %Keep the smallest discrepency so far
        bestSD = sD_cand;
        bestXn = Xn_cand; 
    end
end

fprintf('Best star discrepancy (normalised space) = %.4f\n', bestSD);



% Convert the best point set, calculated above, from the normalised 
% hypercube [0,1]^D into the actual coordinate space used by the problem.
% Each column of 'bestXn' contains values in [0,1], so we apply a
% linear scaling to map them back into their ranges:
%   advertSpend  → [0, 200]
%   targetSpend  → [0, 60]
%   audienceSkew → [0, 1]
% This produces the 20×3 matrix containing evenly spaced sample points

samples = zeros(Nodes,Dimensions); %hold samples
samples(:,1) = minAdvert + bestXn(:,1) * (maxAdvert - minAdvert); %scale
samples(:,2) = minTarget + bestXn(:,2) * (maxTarget - minTarget); %scale
samples(:,3) = minSkew   + bestXn(:,3) * (maxSkew   - minSkew);  %scale

MCSI_samples = exampFunc(samples); 

%===========================================

sliceSkew = 0.5; %Not sure yet but skew which adds the z value of the graph
% which is same length as x and y
%important because "MATLAB fsurf only plots 2 variables, 
%but my function needs 3, so I must create a constant third vector 
%using ones(numel(x),1), multiply by fixedSkew,
%and that gives the 3rd input column."

fsurf(@(x, y) reshape(exampFunc([x(:), y(:),sliceSkew * ones(numel(x),1)]), size(x)), ...
    'MeshDensity', 100, ...
    'ShowContours', 'on', ...
    'LineStyle', ':'); %reshaped and add audience skew because fsurf can only plot for 2 variables
hold on;


xlabel('x');
ylabel('y');
zlabel('exampFunc(x, y)');
title('Surface Plot of exampFunc');


scatter3(samples(:,1),samples(:,2), MCSI_samples, 50, MCSI_samples, 'filled');


hold off;