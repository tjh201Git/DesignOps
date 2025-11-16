%MATLAB BLACKBOX OPTIMIZER
%=============================================================
%Steps
%
% 1) first samples the function using the Star discrepency model
% 2) Plots the result
% 3) Chose a cluster to begin Genetic Algorithm sampling
% 5) repeat until the highest MCSI is found
%====================================================================

minAdvert = 0;
maxAdvert = 200;
minTarget = 0;
maxTarget = 60;
minSkew = 0;
maxSkew = 1;

%===========Star Discrepency Sampling Method===================

Nodes  = 30;
Dimensions  = 3; 
               
% We want to generate many candidate point sets in normalised 3D space
% and measure the star discrepancy of each one. The goal is to find the
% set of points that is the most evenly space, which is the lowest 
% discrepency, to use as our initial sampling design.

numSets = 2000; %Try 2000 different designs
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


%Now plotting the example blackbox function with the 3 variable ranges so
%that we can compare agaisnt the surrogate model 
%its being plotted as a point cloud because its a 4 dimensional plot

Nblackbox = 2000; %Resolution for blackbox plotting

X_Blackbox = [ ... %%
    rand(Nblackbox,1)*(maxAdvert-minAdvert) + minAdvert, ...  % x1: advert
    rand(Nblackbox,1)*(maxTarget-minTarget) + minTarget, ...  % x2: target
    rand(Nblackbox,1)*(maxSkew-minSkew)     + minSkew   ];    % x3: skew

MCSI_true = exampFunc(X_Blackbox);

figure;
scatter3(X_Blackbox(:,1), X_Blackbox(:,2), X_Blackbox(:,3), 25, MCSI_true, 'filled');
hold on;


xlabel('Advert Spending');
ylabel('Target Spending');
zlabel('Audience Skew');
title('4D Plot of the example function');
colorbar;
grid on;
view(45, 25);

%scatter3(samples(:,1),samples(:,2), MCSI_samples, 50, MCSI_samples,'red', 'filled');

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


%========Leave One Out Cross Validation and Quality of Fit===============================
% Leave one out Validation will leave a point out of the surrogate model
% and assess how well the model predicts it and accumulates the squared prediction error.
% This provides an unbiased estimate of surrogate accuracy.

yhatLOO = zeros(Nodes,1); %Initialize leave one out yhat
sseLOO = 0;              %Initialize sum of squared errors  

for i = 1:Nodes 

    %CROSS VALIDATION - use the data of all nodes except the current one
    idx_train = setdiff(1:Nodes, i);  
    
    Phi_train = Phi(idx_train,:);  %trains a phi using the current points except the i'th point
    y_train   = MCSI_samples(idx_train);      
    A_i = pinv(Phi_train) * y_train;
    yhat_i = Phi(i,:) * A_i;   %Makes a surrogate model using every point except i
                               %a different surrogate is made each time in
                               %loop
    yhatLOO(i) = yhat_i;       %Stores the current yhat leave one out value
    
    
    sseLOO = sseLOO + (MCSI_samples(i) - yhat_i)^2;  %calculates the total sum of squared error
                                                     % Accumulate error
                                                     % across all
end

RMSE_LOO = sqrt(sseLOO / Nodes); %Root mean squared of the sum of squared errors

minM = min(MCSI_samples);
maxM = max(MCSI_samples);
rangeM = maxM - minM;
stdM   = std(MCSI_samples);

errorPercent = (RMSE_LOO/rangeM)*100 %the error expressed as a range of the MCSI values

fprintf('LOOCV Error = %.4f (%.2f%% of MCSI range)\n', RMSE_LOO, errorPercent);


%==============Plotting FUll surrogate===================================
figure;
scatter3(Xq(:,1), Xq(:,2), Xq(:,3), 25, yhat, 'filled');
hold on;

% Plot the samples
%scatter3(samples(:,1), samples(:,2), samples(:,3),60, MCSI_samples, 'filled', 'MarkerEdgeColor','k');

xlabel('x1 (Advert)');
ylabel('x2 (Target)');
zlabel('x3 (Skew)');
title('4D Scatter of Surrogate Model (colour = predicted MCSI)');
colorbar;
grid on;
view(45, 25);

%==========Optimizing! Design Space Refinement and Infill==================

origMin = [minAdvert, minTarget, minSkew]; %convienient to have the bounds like this
origMax = [maxAdvert, maxTarget, maxSkew];

%Find the minimum point within the surrogate using fminsearch
% f_sur(x) – Surrogate model evaluation function 
%
% fminsearch requires a function that takes a point x = [x1 x2 x3] in the
% real design space and returns a SINGLE scalar value to minimise.
% This gives a continuous surrogate ŷ(x) that can be optimised using
% fminsearch. The variable 'yhat' from earlier CANNOT be used here because
% yhat only contains surrogate predictions at pre-chosen sample points.
% f_sur(x), in contrast, evaluates the surrogate at ANY point x.
f_sur = @(x) ... 
    [1, ...
     ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1))), ...
     ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2))), ...
     ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3))), ...
     ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1)))^2, ...
     ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2)))^2, ...
     ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3)))^2, ...
     ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1))) * ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2))), ...
     ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1))) * ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3))), ...
     ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2))) * ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3))) ...
    ] * A;
[~, bestIdx] = min(MCSI_samples); %Start the search at the smallest MCSI sample 
x0 = samples(bestIdx,:);           
[xSur_opt, fSur_opt] = fminsearch(f_sur, x0);  %Find the minimum point within the surrogate

xSur_opt = max(xSur_opt, origMin); %Ensure the point isnt outside lower bound
xSur_opt = min(xSur_opt, origMax); %Ensure the point isnt outside upper bound

zoomFactor = 0.5; % ZOOM in by 30% of original size

halfWidth = 0.5 * zoomFactor .* (origMax - origMin); %Finds half the width of the refinement Window

refMin = xSur_opt - halfWidth; %Lower bound of refinement Window
refMax = xSur_opt + halfWidth; %Upper bound of refinement Window

refMin = max(refMin, origMin); %Ensure within bounds
refMax = min(refMax, origMax);

%Refit Surrogate inside the window
%Add 10 new points using the D star in the window and make a surrogate model inside this
%window including the points in this window already from the previous model

%=======NEW surrogate inside the window=====
%Repeat of polynomial surrogate building inside the refinement space

N_infill    = 20;      % number of new points to add
Dimensions  = 3;       
numSets_inf = 2000;     

bestSD_inf  = inf;
bestXn_inf  = [];

for k = 1:numSets_inf
    Xn_cand_inf = rand(N_infill, Dimensions);
    
    sD_cand_inf = starD(Xn_cand_inf);
    
    if sD_cand_inf < bestSD_inf
        bestSD_inf  = sD_cand_inf;
        bestXn_inf  = Xn_cand_inf;
    end
end

%new samples inside the refinement window
samplesInsideWindow = zeros(N_infill, Dimensions);
samplesInsideWindow(:,1) = refMin(1) + bestXn_inf(:,1) * (refMax(1) - refMin(1));  % advert
samplesInsideWindow(:,2) = refMin(2) + bestXn_inf(:,2) * (refMax(2) - refMin(2));  % target
samplesInsideWindow(:,3) = refMin(3) + bestXn_inf(:,3) * (refMax(3) - refMin(3));  % skew

% stick the samples into the blackbox
MCSI_infill = exampFunc(samplesInsideWindow);

samples_new = [samples;samplesInsideWindow]; %Adding new samples to the old
MCSI_samples_new = [MCSI_samples; MCSI_infill]; %Adding new samples to the old
N_new = size(samples_new,1); %new total number of nodes

X_norm_aug = (samples_new - Xmin) ./ (Xmax - Xmin);

x1a = X_norm_aug(:,1);
x2a = X_norm_aug(:,2);
x3a = X_norm_aug(:,3);

Phi_new = [ ...
    ones(N_new,1), ...
    x1a, x2a, x3a, ...
    x1a.^2, x2a.^2, x3a.^2, ...
    x1a.*x2a, x1a.*x3a, x2a.*x3a ];

A_new = pinv(Phi_new) * MCSI_samples_new;  

Xq_new = [ ...
    rand(2000,1)*(maxAdvert-minAdvert) + minAdvert, ...
    rand(2000,1)*(maxTarget-minTarget) + minTarget, ...
    rand(2000,1)*(maxSkew-minSkew)     + minSkew   ];

Xq_norm_new = (Xq_new - Xmin) ./ (Xmax - Xmin); 

x1qn_new = Xq_norm_new(:,1);
x2qn_new = Xq_norm_new(:,2);
x3qn_new = Xq_norm_new(:,3);

Phi_q = [ ...       
    ones(size(Xq_new,1),1), ...
    x1qn_new, x2qn_new, x3qn_new, ...
    x1qn_new.^2, x2qn_new.^2, x3qn_new.^2, ...
    x1qn_new.*x2qn_new, x1qn_new.*x3qn_new, x2qn_new.*x3qn_new ];

yhat_new = Phi_q * A_new;            

%========Plotting new surrogate======================
figure;
scatter3(Xq_new(:,1), Xq_new(:,2), Xq_new(:,3), 25, yhat_new, 'filled');
hold on;



%==========================================================
% Draw refinement window as a transparent 3D box (cube)
%==========================================================
x1_min = refMin(1); x1_max = refMax(1);
x2_min = refMin(2); x2_max = refMax(2);
x3_min = refMin(3); x3_max = refMax(3);

% 8 corner points of the box
C = [ x1_min x2_min x3_min;  % 1
      x1_max x2_min x3_min;  % 2
      x1_max x2_max x3_min;  % 3
      x1_min x2_max x3_min;  % 4
      x1_min x2_min x3_max;  % 5
      x1_max x2_min x3_max;  % 6
      x1_max x2_max x3_max;  % 7
      x1_min x2_max x3_max]; % 8

% Faces of the box (each row indexes into C)
faces = [1 2 3 4;  % bottom (z = min)
         5 6 7 8;  % top    (z = max)
         1 2 6 5;  % side x+
         2 3 7 6;  % side y+
         3 4 8 7;  % side x-
         4 1 5 8]; % side y-

% Plot semi-transparent refinement window
patch('Vertices', C, 'Faces', faces, ...
      'FaceColor', 'none', ...        % no fill, just edges OR:
      'EdgeColor', 'r', ...
      'LineWidth', 1.5);

% If you prefer a light transparent fill, replace FaceColor/FaceAlpha:
% patch('Vertices', C, 'Faces', faces, ...
%       'FaceColor', [1 1 1], 'FaceAlpha', 0.1, ...
%       'EdgeColor', 'k', 'LineWidth', 1.5);
%==========================================================






xlabel('x1 (Advert)');
ylabel('x2 (Target)');
zlabel('x3 (Skew)');
title('4D Scatter of Refined Surrogate Model (colour = predicted MCSI)');
colorbar;
grid on;
view(45, 25);

%========Leave One Out Cross Validation and Quality of Fit of New Surrogate===============================

yhatLOO_new = zeros(N_new,1); %Initialize leave one out yhat
sseLOO_new = 0;              %Initialize sum of squared errors  

for i = 1:N_new 

    %CROSS VALIDATION - use the data of all nodes except the current one
    idx_train = setdiff(1:N_new, i);  
    
    Phi_train = Phi_new(idx_train,:);  %trains a phi using the current points except the i'th point
    y_train   = MCSI_samples_new(idx_train);      
    A_i = pinv(Phi_train) * y_train;
    yhat_i = Phi_new(i,:) * A_i;   %Makes a surrogate model using every point except i
                               %a different surrogate is made each time in
                               %loop
    yhatLOO_new(i) = yhat_i;       %Stores the current yhat leave one out value
    
    
    sseLOO_new = sseLOO_new + (MCSI_samples_new(i) - yhat_i)^2;  %calculates the total sum of squared error
                                                     % Accumulate error
                                                     % across all
end

RMSE_LOO_new = sqrt(sseLOO_new / N_new); %Root mean squared of the sum of squared errors

minM = min(MCSI_samples_new);
maxM = max(MCSI_samples_new);
rangeM = maxM - minM;
stdM   = std(MCSI_samples_new);

errorPercent_new = (RMSE_LOO_new/rangeM)*100 %the error expressed as a range of the MCSI values

fprintf('LOOCV Error For Refined Surrogate = %.4f (%.2f%% of MCSI range)\n', RMSE_LOO_new, errorPercent_new);

%=============Minima of the refined Surrogate====================================
f_sur_new = @(x) ...
    [1, ...
     ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1))), ...
     ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2))), ...
     ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3))), ...
     ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1)))^2, ...
     ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2)))^2, ...
     ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3)))^2, ...
     ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1))) * ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2))), ...
     ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1))) * ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3))), ...
     ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2))) * ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3))) ...
    ] * A_new;

[~, bestIdx] = min(MCSI_samples_new); %Start the search at the smallest new MCSI sample 
x0_new = samples_new(bestIdx,:);         

[MinimumPoint, minimumMCSI] = fminsearch(f_sur_new, x0_new);  % minimise refined surrogate

MinimumPoint = max(MinimumPoint, refMin);
MinimumPoint = min(MinimumPoint, refMax);

fprintf('Refined surrogate minimum at [%.3f, %.3f, %.3f], minimum MCSI = %.4f\n', ...
        MinimumPoint(1), MinimumPoint(2), MinimumPoint(3), minimumMCSI);

%===============Plot the Original MCSI samples VS the NEW=============================

figure;
hold on;

% --- Original samples (all initial star discrepancy points) ---
scatter3(samples(:,1), samples(:,2), samples(:,3), ...
         60, MCSI_samples, 'o', ...
         'MarkerFaceColor', [0 0.447 0.741], ...   % blue
         'MarkerEdgeColor', 'k', ...
         'DisplayName', 'Original Samples');

% --- New infill samples (inside window) ---
scatter3(samplesInsideWindow(:,1), samplesInsideWindow(:,2), samplesInsideWindow(:,3), ...
         80, MCSI_infill, 's', ...
         'MarkerFaceColor', [0.85 0.325 0.098], ... % orange/red
         'MarkerEdgeColor', 'k', ...
         'DisplayName', 'Infill Samples');

xlabel('x1 (Advert)');
ylabel('x2 (Target)');
zlabel('x3 (Skew)');
title('Original vs Infill Sample Points (coloured by MCSI)');
colorbar;
grid on;
view(45,25);
legend('Location','best');
hold off;
%=======RBF Surrogate=====