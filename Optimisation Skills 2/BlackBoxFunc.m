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

%===========Star Discrepency Sampling Method===================

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

%=================Plotting==========================

%sliceSkew = 0.5; %Not sure yet but skew which adds the z value of the graph
%sliceSkew * ones(numel(x),1)

% which is same length as x and y
%important because "MATLAB fsurf only plots 2 variables, 
%but my function needs 3, so I must create a constant third vector 
%using ones(numel(x),1), multiply by fixedSkew,
%and that gives the 3rd input column."

fsurf(@(x, y) reshape(exampFunc([x(:), y(:)]), size(x)), ...
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


%===========Build a Surrogate Model for Fitting==============================


%=======Polynomial Surrogate=====

%Second order Polynomial in 3 variables
%   ŷ = a0
%      + a1*x1 + a2*x2 + a3*x3      (linear terms)
%      + a4*x1^2 + a5*x2^2 + a6*x3^2   (quadratic terms)
%      + a7*x1*x2 + a8*x1*x3 + a9*x2*x3  (interaction terms)
% 
% 10 unknown coefficients (a0 ... a9)
% 
% Determine the unknown coefficients with linear regression
%
%       y ≈ Φ * A
%
% where:
%   Φ is the matrix of basis functions,
%   A is the coefficient vector,
%   y is the vector of MCSI outputs.
%

%Fit the surrogate on normalized inputs
Xmin = [minAdvert, minTarget, minSkew]; % Normalise samples to [0,1]
Xmax = [maxAdvert, maxTarget, maxSkew];

X_norm = (samples - Xmin) ./ (Xmax - Xmin);  % 20x3 matrix in [0,1]​

x1 = X_norm(:,1); %extract each input dimension from the 20
x2 = X_norm(:,2);
x3 = X_norm(:,3);

Phi = [ ... % Phi = [1, x1, x2, x3, x1^2, x2^2, x3^2, x1*x2, x1*x3, x2*x3]
    ones(Nodes,1), ...   
    x1, x2, x3, ...               
    x1.^2, x2.^2, x3.^2, ...       
    x1.*x2, x1.*x3, x2.*x3];       

A = pinv(Phi) * MCSI_samples; % Solve for coefficients using pseudo-inverse
                              % find the least squares solution that fits

%Build query points and visualise the estimated predicted polynomial                              

%First Create a 2000 by 3 matrix to visualise the points(2000 is a
%resolution)
Xq = [ ...
    rand(2000,1)*(maxAdvert-minAdvert) + minAdvert, ...
    rand(2000,1)*(maxTarget-minTarget) + minTarget, ...
    rand(2000,1)*(maxSkew-minSkew)     + minSkew   ];


Xq_norm = (Xq - Xmin) ./ (Xmax - Xmin); % Normalise the query points to 0-1
%Surrogate modelling in a normalised space improves fit quality

x1qn = Xq_norm(:,1);
x2qn = Xq_norm(:,2);
x3qn = Xq_norm(:,3);

Phi_q = [ ...        %Constrauct the same Phi matrix for the query points
    ones(size(Xq,1),1), ...
    x1qn, x2qn, x3qn, ...
    x1qn.^2, x2qn.^2, x3qn.^2, ...
    x1qn.*x2qn, x1qn.*x3qn, x2qn.*x3qn ];

%Evaluate the surrogate model at each query pont 
yhat = Phi_q * A;                     % nq^2 x 1

figure;
scatter3(Xq(:,1), Xq(:,2), Xq(:,3), 25, yhat, 'filled');
hold on;

% Plot the samples
scatter3(samples(:,1), samples(:,2), samples(:,3),60, MCSI_samples, 'filled', 'MarkerEdgeColor','k');

xlabel('x1 (Advert)');
ylabel('x2 (Target)');
zlabel('x3 (Skew)');
title('4D Scatter of Surrogate Model (colour = predicted MCSI)');
colorbar;
grid on;
view(45, 25);

%=======RBF Surrogate=====