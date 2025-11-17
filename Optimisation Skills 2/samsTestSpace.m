% sam test space
clear all

minAdvert = 0;
maxAdvert = 200;
minTarget = 0;
maxTarget = 60;
minSkew = 0;
maxSkew = 1;
maxVariables = [maxAdvert, maxTarget, maxSkew];
minVariables = [minAdvert, minTarget, minSkew];


dimensions = 3;
iterations = 1000;
nodes = 30;
saveDStarDistribution = false;

blackBoxFunc = @griewank3;

[pointsMatrix, bestDStar] = bestDStarPointDistribution(dimensions, nodes, iterations, saveDStarDistribution);



% denormalisedHyperCube = zeros(Nodes,Dimensions); %hold samples
% denormalisedHyperCube(:,1) = denormaliseMatrix(pointsMatrix(:,1), minAdvert, maxAdvert); %scale
% denormalisedHyperCube(:,2) = denormaliseMatrix(pointsMatrix(:,2), minTarget, maxTarget); %scale
% denormalisedHyperCube(:,3) = denormaliseMatrix(pointsMatrix(:,3), minSkew, maxSkew); %scale

denormalisedHyperCube = denormaliseMatrix(pointsMatrix, minVariables, maxVariables);



%================= Interactive sample taking 1st round ==========================
MCSI_samples = blackBoxFunc(denormalisedHyperCube); % or take interactive samples
%================================================================================

% plot the first samples
figure('Name', 'Plot of the First Run', 'NumberTitle', 'off');
scatter3(denormalisedHyperCube(:,1), denormalisedHyperCube(:,2), denormalisedHyperCube(:,3), 25, MCSI_samples, 'filled');
hold on;

xlabel('Advert Spending');
ylabel('Target Spending');
zlabel('Audience Skew');
title('4D Plot of the example function');
colorbar;
grid on;
view(45, 25);
hold off;



%===========Build a Surrogate Model for Fitting==============================

normalisedSamples = normaliseMatrix(denormalisedHyperCube, minVariables, maxVariables);

%extract each input dimension
advertWeights = normalisedSamples(:,1); 
targetWeights = normalisedSamples(:,2);
skewWeights = normalisedSamples(:,3);

%==============Decide on the Polynomial Order==============================

polynomial_orders_to_test = 5;

phi_orders = zeros(polynomial_orders_to_test, 1);

for i = 1:polynomial_orders_to_test
    phi_orders(i) = buildPhi(nodes,i,advertWeights,targetWeights,skewWeights);
end

%========Leave One Out Cross Validation and Quality of Fit===============================
% Leave one out Validation will leave a point out of the surrogate model
% and assess how well the model predicts it and accumulates the squared prediction error.
% This provides an unbiased estimate of surrogate accuracy.

RMSE_orders = zeros(polynomial_orders_to_test,1);

for i = 1:polynomial_orders_to_test
    RMSE_orders(i) = LOOV(nodes,phi_orders(i),MCSI_samples);
end

rangeM = max(MCSI_samples) - min(MCSI_samples);

[bestRMSE, besPolynomialOrder] = min(RMSE_orders);

errorPercent = (bestRMSE/rangeM)*100; %the error expressed as a range of the MCSI values

fprintf('POLY First Surrogate LOOCV Error = %.4f (%.2f%% of MCSI range), Order:%.1f\n', bestRMSE, errorPercent, BestOrder);

phi = phi_orders(besPolynomialOrder);

A = pinv(phi) * MCSI_samples; % Solve for coefficients using pseudo-inverse
                              % find the least squares solution that fits


%Build query points and visualise the estimated predicted polynomial                              
%First Create a 2000 by 3 matrix to visualise the points(2000 is a
%resolution)
% Xq = [ ...
%     rand(2000,1)*(maxAdvert-minAdvert) + minAdvert, ...
%     rand(2000,1)*(maxTarget-minTarget) + minTarget, ...
%     rand(2000,1)*(maxSkew-minSkew)     + minSkew   ];

queryPointsNum = 2000;
normalisedQueryPoints = [rand(queryPointsNum, 1); rand(queryPointsNum, 1); rand(queryPointsNum, 1)];
denormalisedQueryPoints = denormaliseMatrix(normalisedQueryPoints, minVariables, maxVariables);

advertWeightsQuery = normalisedQueryPoints(:,1); 
targetWeightsQuery = normalisedQueryPoints(:,2);
skewWeightsQuery = normalisedQueryPoints(:,3);

phiQuery = buildPhi_q(denormalisedQueryPoints, besPolynomialOrder, advertWeightsQuery, targetWeightsQuery, skewWeightsQuery);

%Evaluate the surrogate model at each query pont 
yhat = phiQuery * A;                     % nq^2 x 1

% %==============Plotting FUll surrogate===================================
figure('Name', 'Plot of the Surrogate Model', 'NumberTitle', 'off');

scatter3(advertWeightsQuery, targetWeightsQuery, skewWeightsQuery, 25, yhat, 'filled');
hold on;

% Plot the samples
%scatter3(samples(:,1), samples(:,2), samples(:,3),60, MCSI_samples, 'filled', 'MarkerEdgeColor','k');

xlabel('x1 (Advert)');
ylabel('x2 (Target)');
zlabel('x3 (Skew)');
title('Poly 4D Scatter of Surrogate Model (colour = predicted MCSI)');
colorbar;
grid on;
view(45, 25);

% %==========Optimizing! Design Space Refinement and Infill=====================
% 
% origMin = [minAdvert, minTarget, minSkew]; %convienient to have the bounds like this
% origMax = [maxAdvert, maxTarget, maxSkew];
% 
% %Find the minimum point within the surrogate using Genetic Algorithm
% % GA is used here because it handles bound constraints and
% % performs a global search which reduces the chance of becoming trapped in a
% % local minimum of the surrogate. Unlike fminsearch (which is
% % unconstrained), GA ensures that all candidate points remain inside the
% % original design space
% 
% f_sur = @(x) ...  %Builds  callable version of the polynomial surrogate to
%     [1, ...       %to be put into the Genetic Algorithm function. It takes
%      ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1))), ...%any point in the space and 
%      ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2))), ...%normalises it to 0-1 and 
%      ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3))), ... %reconstructs [1, x1, x2, x3, x1^2, x2^2, x3^2, x1*x2, x1*x3, x2*x3]
%      ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1)))^2, ...
%      ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2)))^2, ...
%      ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3)))^2, ...
%      ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1))) * ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2))), ...
%      ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1))) * ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3))), ...
%      ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2))) * ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3))) ...
%     ] * A;
% 
% lowerBound = Xmin;
% upperBound = Xmax;
% 
% [xSur_opt, fSur_opt] = ga(@(x) f_sur(x), 3, [], [], [], [], lowerBound, upperBound, [], optsGA);
% 
% xSur_opt = max(xSur_opt, origMin); %Ensure the point isnt outside lower bound
% xSur_opt = min(xSur_opt, origMax); %Ensure the point isnt outside upper bound
% 
% zoomFactor = 0.5; % ZOOM in by 30% of original size
% 
% halfWidth = 0.5 * zoomFactor .* (origMax - origMin); %Finds half the width of the refinement Window
% 
% refMin = xSur_opt - halfWidth; %Lower bound of refinement Window
% refMax = xSur_opt + halfWidth; %Upper bound of refinement Window
% 
% refMin = max(refMin, origMin); %Ensure within bounds
% refMax = min(refMax, origMax);
% 
% %Refit Surrogate inside the window
% %Add 10 new points using the D star in the window and make a surrogate model inside this
% %window including the points in this window already from the previous model
% 
% %==============NEW surrogate inside the window===========
% %Repeat of polynomial surrogate building inside the refinement space
% 
% N_infill    = PolysurrogateNodes;      % number of new points to add
% Dimensions  = 3;       
% numSets_inf = 2000;     
% 
% bestSD_inf  = inf;
% bestXn_inf  = [];
% 
% for k = 1:numSets_inf
%     Xn_cand_inf = rand(N_infill, Dimensions);
% 
%     sD_cand_inf = starD(Xn_cand_inf);
% 
%     if sD_cand_inf < bestSD_inf
%         bestSD_inf  = sD_cand_inf;
%         bestXn_inf  = Xn_cand_inf;
%     end
% end
% 
% %new samples inside the refinement window
% samplesInsideWindow = zeros(N_infill, Dimensions);
% samplesInsideWindow(:,1) = refMin(1) + bestXn_inf(:,1) * (refMax(1) - refMin(1));  % advert
% samplesInsideWindow(:,2) = refMin(2) + bestXn_inf(:,2) * (refMax(2) - refMin(2));  % target
% samplesInsideWindow(:,3) = refMin(3) + bestXn_inf(:,3) * (refMax(3) - refMin(3));  % skew
% 
% % stick the samples into the blackbox
% MCSI_infill = ChosenFunc(samplesInsideWindow);
% 
% samples_new = [samples;samplesInsideWindow]; %Adding new samples to the old
% MCSI_samples_new = [MCSI_samples; MCSI_infill]; %Adding new samples to the old
% N_new = size(samples_new,1); %new total number of nodes
% 
% X_norm_aug = (samples_new - Xmin) ./ (Xmax - Xmin);
% 
% x1a = X_norm_aug(:,1);
% x2a = X_norm_aug(:,2);
% x3a = X_norm_aug(:,3);
% 
% currentNodeCount = PolysurrogateNodes + Nodes;
% 
% %===========Find the best Order=========================
% Phi1order = buildPhi(currentNodeCount,1,x1a,x2a,x3a);
% Phi2order = buildPhi(currentNodeCount,2,x1a,x2a,x3a);
% Phi3order = buildPhi(currentNodeCount,3,x1a,x2a,x3a);
% Phi4order = buildPhi(currentNodeCount,4,x1a,x2a,x3a);
% Phi5order = buildPhi(currentNodeCount,5,x1a,x2a,x3a);
% 
% RMSE_orders = zeros(5,1);
% 
% RMSE_orders(1) = LOOV(currentNodeCount,Phi1order,MCSI_samples_new); %Root mean squared of the sum of squared errors
% RMSE_orders(2) = LOOV(currentNodeCount,Phi2order,MCSI_samples_new); %Root mean squared of the sum of squared errors
% RMSE_orders(3) = LOOV(currentNodeCount,Phi3order,MCSI_samples_new); %Root mean squared of the sum of squared errors
% RMSE_orders(4) = LOOV(currentNodeCount,Phi4order,MCSI_samples_new); %Root mean squared of the sum of squared errors
% RMSE_orders(5) = LOOV(currentNodeCount,Phi5order,MCSI_samples_new); %Root mean squared of the sum of squared errors
% 
% rangeM = max(MCSI_samples_new) - min(MCSI_samples_new);
% 
% [bestRMSE, BestOrder] = min(RMSE_orders);
% 
% errorPercent = (bestRMSE/rangeM)*100; %the error expressed as a range of the MCSI values
% 
% fprintf('POLY Second Surrogate LOOCV Error = %.4f (%.2f%% of MCSI range), Order:%.1f\n', bestRMSE, errorPercent, BestOrder);
% 
% switch BestOrder %Switch case for which best order was chosen
%     case 1
%         Phi_new = Phi1order;
%     case 2
%         Phi_new = Phi2order;
%     case 3
%         Phi_new = Phi3order;
%     case 4
%         Phi_new = Phi4order;
%     case 5
%         Phi_new = Phi5order;
%     otherwise
%         error('BestOrder must be an integer between 1 and 5.');
% end
% 
% A_new = pinv(Phi_new) * MCSI_samples_new;  
% 
% 
% Xq_new = [ ...
%     rand(2000,1)*(maxAdvert-minAdvert) + minAdvert, ...
%     rand(2000,1)*(maxTarget-minTarget) + minTarget, ...
%     rand(2000,1)*(maxSkew-minSkew)     + minSkew   ];
% 
% Xq_norm_new = (Xq_new - Xmin) ./ (Xmax - Xmin); 
% 
% x1qn_new = Xq_norm_new(:,1);
% x2qn_new = Xq_norm_new(:,2);
% x3qn_new = Xq_norm_new(:,3);
% 
% Phi1order_q = buildPhi_q(Xq_new,1,x1qn_new,x2qn_new,x3qn_new);
% Phi2order_q = buildPhi_q(Xq_new,2,x1qn_new,x2qn_new,x3qn_new);
% Phi3order_q = buildPhi_q(Xq_new,3,x1qn_new,x2qn_new,x3qn_new);
% Phi4order_q = buildPhi_q(Xq_new,4,x1qn_new,x2qn_new,x3qn_new);
% Phi5order_q = buildPhi_q(Xq_new,5,x1qn_new,x2qn_new,x3qn_new);
% 
% switch BestOrder %Switch case for which best order was chosen
%     case 1
%         Phi_q_new = Phi1order_q;
%     case 2
%         Phi_q_new = Phi2order_q;
%     case 3
%         Phi_q_new = Phi3order_q;
%     case 4
%         Phi_q_new = Phi4order_q;
%     case 5
%         Phi_q_new = Phi5order_q;
%     otherwise
%         error('BestOrder must be an integer between 1 and 5.');
% end
% 
% yhat_new = Phi_q_new * A_new;            
% 
% %========Plotting new surrogate======================
% figure('Name', 'Plot of the New Surrogate (with Refinement)', 'NumberTitle', 'off');
% 
% scatter3(Xq_new(:,1), Xq_new(:,2), Xq_new(:,3), 25, yhat_new, 'filled');
% hold on;
% 
% %====================Draw refinement window=====================
% 
% x1_min = refMin(1); x1_max = refMax(1);
% x2_min = refMin(2); x2_max = refMax(2);
% x3_min = refMin(3); x3_max = refMax(3);
% 
% % 8 corner points of the box
% C = [ x1_min x2_min x3_min;  % 1
%       x1_max x2_min x3_min;  % 2
%       x1_max x2_max x3_min;  % 3
%       x1_min x2_max x3_min;  % 4
%       x1_min x2_min x3_max;  % 5
%       x1_max x2_min x3_max;  % 6
%       x1_max x2_max x3_max;  % 7
%       x1_min x2_max x3_max]; % 8
% 
% % Faces of the box (each row indexes into C)
% faces = [1 2 3 4;  % bottom (z = min)
%          5 6 7 8;  % top    (z = max)
%          1 2 6 5;  % side x+
%          2 3 7 6;  % side y+
%          3 4 8 7;  % side x-
%          4 1 5 8]; % side y-
% 
% % Plot semi-transparent refinement window
% patch('Vertices', C, 'Faces', faces, ...
%       'FaceColor', 'none', ...        % no fill, just edges OR:
%       'EdgeColor', 'r', ...
%       'LineWidth', 1.5);
% 
% %==========================================================
% 
% xlabel('x1 (Advert)');
% ylabel('x2 (Target)');
% zlabel('x3 (Skew)');
% title('POLY 4D Scatter of Refined Surrogate Model (colour = predicted MCSI)');
% colorbar;
% grid on;
% view(45, 25);
% 
% %=============Minima of the refined Surrogate====================================
% 
% %=====Run a Genetic Algorithim on the Surrogate to find Minima========
% 
% f_sur_new = @(x) ...
%     [1, ...
%      ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1))), ...
%      ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2))), ...
%      ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3))), ...
%      ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1)))^2, ...
%      ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2)))^2, ...
%      ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3)))^2, ...
%      ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1))) * ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2))), ...
%      ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1))) * ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3))), ...
%      ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2))) * ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3))) ...
%     ] * A_new;
% 
% lowerBound = Xmin;
% upperBound = Xmax;
% 
% [xGA, fGA] = ga(@(x) f_sur_new(x), 3, [], [], [], [], lowerBound, upperBound, [], optsGA);
% 
% %trueAtGA = ChosenFunc(xGA); 
% 
% 
% fprintf('GA surrogate min at [%.3f, %.3f, %.3f]\n', xGA(1), xGA(2), xGA(3));
% %fprintf('  Surrogate MCSI  = %.4f\n', fGA);
% 
% 
% %===============Third Surrogate=============
% %Will now use the GA optimum of the Second surrogate as the center of the
% %smaller refinement window. This Repeats the design refinement logic but
% %with a tighter zoom
% 
% zoomFactor3 = 0.2;  
% 
% center3 = max(min(xGA, origMax), origMin);
% 
% halfWidth3 = 0.5 * zoomFactor3 .* (origMax - origMin);
% 
% refMin3 = center3 - halfWidth3;
% refMax3 = center3 + halfWidth3;
% 
% refMin3 = max(refMin3, origMin);
% refMax3 = min(refMax3, origMax);
% 
% N_infill3    = PolythirdSurrogateNodes;    % number of new points to add (3rd refinement)
% numSets_inf3 = 2000;
% 
% bestSD_inf3  = inf;
% bestXn_inf3  = [];
% 
% for k = 1:numSets_inf3
%     Xn_cand_inf3 = rand(N_infill3, Dimensions);   
%     sD_cand_inf3 = starD(Xn_cand_inf3);          
% 
%     if sD_cand_inf3 < bestSD_inf3
%         bestSD_inf3  = sD_cand_inf3;
%         bestXn_inf3  = Xn_cand_inf3;
%     end
% end
% 
% samplesInsideWindow3 = zeros(N_infill3, Dimensions);
% samplesInsideWindow3(:,1) = refMin3(1) + bestXn_inf3(:,1) * (refMax3(1) - refMin3(1));  % advert
% samplesInsideWindow3(:,2) = refMin3(2) + bestXn_inf3(:,2) * (refMax3(2) - refMin3(2));  % target
% samplesInsideWindow3(:,3) = refMin3(3) + bestXn_inf3(:,3) * (refMax3(3) - refMin3(3));  % skew
% 
% MCSI_infill3 = ChosenFunc(samplesInsideWindow3);% Evaluate black-box on the 3rd refinement points
% 
% % add all data: original + 2nd refinement + 3rd refinement
% samples_third      = [samples_new; samplesInsideWindow3]; 
% MCSI_samples_third = [MCSI_samples_new; MCSI_infill3];
% N_third            = size(samples_third, 1);
% 
% X_norm_third = (samples_third - Xmin) ./ (Xmax - Xmin);
% 
% x1t = X_norm_third(:,1);
% x2t = X_norm_third(:,2);
% x3t = X_norm_third(:,3);
% 
% Phi1order = buildPhi(N_third,1,x1t,x2t,x3t);
% Phi2order = buildPhi(N_third,2,x1t,x2t,x3t);
% Phi3order = buildPhi(N_third,3,x1t,x2t,x3t);
% Phi4order = buildPhi(N_third,4,x1t,x2t,x3t);
% Phi5order = buildPhi(N_third,5,x1t,x2t,x3t);
% 
% RMSE_orders = zeros(5,1);
% 
% RMSE_orders(1) = LOOV(N_third,Phi1order,MCSI_samples_third); %Root mean squared of the sum of squared errors
% RMSE_orders(2) = LOOV(N_third,Phi2order,MCSI_samples_third); %Root mean squared of the sum of squared errors
% RMSE_orders(3) = LOOV(N_third,Phi3order,MCSI_samples_third); %Root mean squared of the sum of squared errors
% RMSE_orders(4) = LOOV(N_third,Phi4order,MCSI_samples_third); %Root mean squared of the sum of squared errors
% RMSE_orders(5) = LOOV(N_third,Phi5order,MCSI_samples_third); %Root mean squared of the sum of squared errors
% 
% rangeM = max(MCSI_samples_third) - min(MCSI_samples_third);
% 
% [bestRMSE, BestOrder] = min(RMSE_orders);
% 
% errorPercent = (bestRMSE/rangeM)*100; %the error expressed as a range of the MCSI values
% 
% fprintf('POLY Secnd Surrogate LOOCV Error = %.4f (%.2f%% of MCSI range), Order:%.1f\n', bestRMSE, errorPercent, BestOrder);
% 
% switch BestOrder %Switch case for which best order was chosen
%     case 1
%         Phi_third = Phi1order;
%     case 2
%         Phi_third = Phi2order;
%     case 3
%         Phi_third = Phi3order;
%     case 4
%         Phi_third = Phi4order;
%     case 5
%         Phi_third = Phi5order;
%     otherwise
%         error('BestOrder must be an integer between 1 and 5.');
% end
% 
% A_third = pinv(Phi_third) * MCSI_samples_third;
% 
% Xq_third = [ ...
%     rand(2000,1)*(maxAdvert-minAdvert) + minAdvert, ...
%     rand(2000,1)*(maxTarget-minTarget) + minTarget, ...
%     rand(2000,1)*(maxSkew-minSkew)     + minSkew   ];
% 
% Xq_norm_third = (Xq_third - Xmin) ./ (Xmax - Xmin);
% 
% x1q3 = Xq_norm_third(:,1);
% x2q3 = Xq_norm_third(:,2);
% x3q3 = Xq_norm_third(:,3);
% 
% Phi1order_q = buildPhi_q(Xq_third,1,x1q3,x2q3,x3q3);
% Phi2order_q = buildPhi_q(Xq_third,2,x1q3,x2q3,x3q3);
% Phi3order_q = buildPhi_q(Xq_third,3,x1q3,x2q3,x3q3);
% Phi4order_q = buildPhi_q(Xq_third,4,x1q3,x2q3,x3q3);
% Phi5order_q = buildPhi_q(Xq_third,5,x1q3,x2q3,x3q3);
% 
% switch BestOrder %Switch case for which best order was chosen
%     case 1
%         Phi_q3 = Phi1order_q;
%     case 2
%         Phi_q3 = Phi2order_q;
%     case 3
%         Phi_q3 = Phi3order_q;
%     case 4
%         Phi_q3 = Phi4order_q;
%     case 5
%         Phi_q3 = Phi5order_q;
%     otherwise
%         error('BestOrder must be an integer between 1 and 5.');
% end
% 
% yhat_third = Phi_q3 * A_third;
% 
% figure('Name', 'Third Refined Surrogate (zoom = 0.3)', 'NumberTitle', 'off');
% scatter3(Xq_third(:,1), Xq_third(:,2), Xq_third(:,3), 25, yhat_third, 'filled');
% hold on;
% 
% %Quality of Fit for the Third Surrogate
% yhatLOO_third = zeros(N_third,1);
% sseLOO_third  = 0;
% 
% for i = 1:N_third
%     idx_train = setdiff(1:N_third, i);
%     Phi_train = Phi_third(idx_train,:);
%     y_train   = MCSI_samples_third(idx_train);
% 
%     A_i = pinv(Phi_train) * y_train;
%     yhat_i = Phi_third(i,:) * A_i;
% 
%     yhatLOO_third(i) = yhat_i;
%     sseLOO_third = sseLOO_third + (MCSI_samples_third(i) - yhat_i)^2;
% end
% 
% RMSE_LOO_third = sqrt(sseLOO_third / N_third);
% minM3 = min(MCSI_samples_third);
% maxM3 = max(MCSI_samples_third);
% rangeM3 = maxM3 - minM3;
% 
% errorPercent_third = (RMSE_LOO_third / rangeM3)*100;
% 
% fprintf('LOOCV Error For THIRD Refined Surrogate = %.4f (%.2f%% of MCSI range)\n', RMSE_LOO_third, errorPercent_third);
% 
% %Genetic Algorithim for the Third Surrogate
% f_sur_third = @(x) ...
%     [1, ...
%      ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1))), ...
%      ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2))), ...
%      ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3))), ...
%      ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1)))^2, ...
%      ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2)))^2, ...
%      ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3)))^2, ...
%      ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1))) * ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2))), ...
%      ((x(1)-Xmin(1))/(Xmax(1)-Xmin(1))) * ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3))), ...
%      ((x(2)-Xmin(2))/(Xmax(2)-Xmin(2))) * ((x(3)-Xmin(3))/(Xmax(3)-Xmin(3))) ...
%     ] * A_third;
% 
% [xGA_third, fGA_third] = ga(@(x) f_sur_third(x), 3, [], [], [], [], ...
%                             lowerBound, upperBound, [], optsGA);
% 
% %trueAtGA_third = ChosenFunc(xGA_third);
% %fprintf('  True MCSI at that point = %.4f\n', trueAtGA_third);
% 
% fprintf('3rd Surrogate GA min at [%.3f, %.3f, %.3f]\n', xGA_third(1), xGA_third(2), xGA_third(3));
% fprintf('  Surrogate MCSI (3rd) = %.4f\n', fGA_third);
% 
% %==============Plotting final surrogate with refined window======================
% figure('Name', 'Plot of the New Surrogate (with Refinement)', 'NumberTitle', 'off');
% 
% scatter3(Xq_third(:,1), Xq_third(:,2), Xq_third(:,3), 25, yhat_third, 'filled');
% hold on;
% 
% %====================Draw refinement window=====================
% 
% x1_min = refMin(1); x1_max = refMax(1);
% x2_min = refMin(2); x2_max = refMax(2);
% x3_min = refMin(3); x3_max = refMax(3);
% 
% % 8 corner points of the box
% C = [ x1_min x2_min x3_min;  % 1
%       x1_max x2_min x3_min;  % 2
%       x1_max x2_max x3_min;  % 3
%       x1_min x2_max x3_min;  % 4
%       x1_min x2_min x3_max;  % 5
%       x1_max x2_min x3_max;  % 6
%       x1_max x2_max x3_max;  % 7
%       x1_min x2_max x3_max]; % 8
% 
% % Faces of the box (each row indexes into C)
% faces = [1 2 3 4;  % bottom (z = min)
%          5 6 7 8;  % top    (z = max)
%          1 2 6 5;  % side x+
%          2 3 7 6;  % side y+
%          3 4 8 7;  % side x-
%          4 1 5 8]; % side y-
% 
% % Plot semi-transparent refinement window
% patch('Vertices', C, 'Faces', faces, ...
%       'FaceColor', 'none', ...        % no fill, just edges OR:
%       'EdgeColor', 'r', ...
%       'LineWidth', 1.5);
% 
% x1_min3 = refMin3(1); x1_max3 = refMax3(1);
% x2_min3 = refMin3(2); x2_max3 = refMax3(2);
% x3_min3 = refMin3(3); x3_max3 = refMax3(3);
% 
% % 8 corner points of the box
% C = [ x1_min3 x2_min3 x3_min3;  % 1
%       x1_max3 x2_min3 x3_min3;  % 2
%       x1_max3 x2_max3 x3_min3;  % 3
%       x1_min3 x2_max3 x3_min3;  % 4
%       x1_min3 x2_min3 x3_max3;  % 5
%       x1_max3 x2_min3 x3_max3;  % 6
%       x1_max3 x2_max3 x3_max3;  % 7
%       x1_min3 x2_max3 x3_max3]; % 8
% 
% % Faces of the box (each row indexes into C)
% faces3 = [1 2 3 4;  % bottom (z = min)
%          5 6 7 8;  % top    (z = max)
%          1 2 6 5;  % side x+
%          2 3 7 6;  % side y+
%          3 4 8 7;  % side x-
%          4 1 5 8]; % side y-
% 
% % Plot semi-transparent refinement window
% patch('Vertices', C, 'Faces', faces3, ...
%       'FaceColor', 'none', ...        % no fill, just edges OR:
%       'EdgeColor', 'g', ...
%       'LineWidth', 1.5);
% 
% %==========================================================
% 
% xlabel('x1 (Advert)');
% ylabel('x2 (Target)');
% zlabel('x3 (Skew)');
% title('POLY 4D Scatter of Refined Surrogate Model (colour = predicted MCSI)');
% colorbar;
% grid on;
% view(45, 25);
% 
% %===============Plot the Original MCSI samples VS the NEW=============================
% figure('Name', 'Original Samples Vs Refined1 and Refined2', 'NumberTitle', 'off');
% hold on;
% 
% % --- Original samples (all initial star discrepancy points) ---
% scatter3(samples(:,1), samples(:,2), samples(:,3), ...
%          60, MCSI_samples, 'o', ...
%          'MarkerFaceColor', [0 0.447 0.741], ...   % blue
%          'MarkerEdgeColor', 'k', ...
%          'DisplayName', 'Original Samples');
% 
% % --- New infill samples (inside window) ---
% scatter3(samplesInsideWindow(:,1), samplesInsideWindow(:,2), samplesInsideWindow(:,3), ...
%          80, MCSI_infill, 's', ...
%          'MarkerFaceColor', [0.85 0.325 0.098], ... % orange/red
%          'MarkerEdgeColor', 'k', ...
%          'DisplayName', 'Infill Samples');
% 
% % --- New infill samples (inside window) ---
% scatter3(samplesInsideWindow3(:,1), samplesInsideWindow3(:,2), samplesInsideWindow3(:,3), ...
%          80, MCSI_infill, 's', ...
%          'MarkerFaceColor', [0.4660 0.6740 0.1880], ... % orange/red
%          'MarkerEdgeColor', 'k', ...
%          'DisplayName', 'Infill Samples');
% 
% xlabel('x1 (Advert)');
% ylabel('x2 (Target)');
% zlabel('x3 (Skew)');
% title('POLY Original vs Infill Sample Points (coloured by MCSI)');
% colorbar;
% grid on;
% view(45,25);
% legend('Location','best');
% hold off;
% 
% 
% %=================RBF Surrogate==========================================
% %Now building an RBF Surrogate using a Gaussian Basis
% %   psi(r) = exp(-(r/L)^2)
% % where r is Euclidean distance in *normalised* space, and L is the
% % length scale (shape parameter).
% 
% %X_norm;                % Defined when making the Polynomial Surrogate
% %Y  = MCSI_samples;          % Defined when making the Polynomial Surrogate
% 
% 
% N  = size(X_norm,1);            
% 
% %Choose a length scale L (Franke's formula, in normalised space)
% %    L = D / (0.8 * sqrt(N)^(1/4))
% %    where D is the "diameter" of the sample cloud.
% D = norm(max(X_norm) - min(X_norm));                 % approximate diameter in [0,1]^3
% L = D / (0.8 * (sqrt(N))^(1/4));             % Equation from the notes
% 
% 
% % Build the NxN Phi matrix of pairwise Gaussian RBF values
% % phi(i,j) = psi( || X_norm(i,:) - X_norm(j,:) || )
% Phi_RBF = zeros(N,N);
% for i = 1:N
%     for j = 1:N
%         r_ij = sqrt( sum( (X_norm(i,:) - X_norm(j,:)).^2 , 2) );  % Euclidean distance
%         Phi_RBF(i,j) = exp( - (r_ij / L)^2 );             % Gaussian RBF
%     end
% end
% 
% % Solve for weights W in:  Y = Phi_RBF * W
% %(interpolates sample points exactly)
% W_RBF = Phi_RBF \ MCSI_samples;    % same as inv(Phi_RBF)*Y but more stable
% 
% Nq = size(Xq_norm,1); %Using Xq and Xq_Norm from the previous bit of Code
% yhat_RBF = zeros(Nq,1);
% 
% for i = 1:Nq
%     % distance from query point i to all sample points in normalised space
%     r_i = sqrt( sum( (X_norm - Xq_norm(i,:)).^2 , 2) );   % N×1
%     psi_i = exp( - (r_i / L).^2 );                    % N×1
%     % RBF prediction is weighted sum of these Gaussians
%     yhat_RBF(i) = psi_i' * W_RBF;                     % scalar
% end
% 
% figure('Name','RBF Surrogate Model','NumberTitle','off');
% scatter3(Xq(:,1), Xq(:,2), Xq(:,3), 25, yhat_RBF, 'filled');
% hold on;
% 
% % Optionally overlay the actual sample points:
% scatter3(samples(:,1), samples(:,2), samples(:,3), ...
%          60, MCSI_samples, 'o', 'MarkerEdgeColor','k', ...
%          'MarkerFaceColor','w');
% 
% xlabel('x1 (Advert)');
% ylabel('x2 (Target)');
% zlabel('x3 (Skew)');
% title('RBF 4D Scatter of RBF Surrogate (colour = predicted MCSI)');
% colorbar;
% grid on;
% view(45,25);
% hold off;
% 
% %==================Optimize!================================================
% 
% %First Optimize Zoom by factor and Infill
% 
% f_RBF = @(x) ...
%     exp( - ( sqrt( sum( (X_norm - ...
%               ((x - Xmin) ./ (Xmax - Xmin)) ).^2 , 2) ) / L ).^2 )' ...
%     * W_RBF;
% 
% [xRBF_opt, fRBF_opt] = ga(@(x) f_RBF(x), 3, [], [], [], [], lowerBound, upperBound, [], optsGA);
% 
% xRBF_opt = max(xRBF_opt, origMin); %ensure its inside the bounds
% xRBF_opt = min(xRBF_opt, origMax);
% 
% fprintf('RBF surrogate GA min at [%.3f, %.3f, %.3f], surrogate = %.4f\n',xRBF_opt(1), xRBF_opt(2), xRBF_opt(3), fRBF_opt);
% 
% zoomFactorRBF = 0.5;  % 50% of original box size
% 
% halfWidthRBF = 0.5 * zoomFactorRBF .* (origMax - origMin);
% 
% refMinRBF = xRBF_opt - halfWidthRBF;
% refMaxRBF = xRBF_opt + halfWidthRBF;
% 
% refMinRBF = max(refMinRBF, origMin);
% refMaxRBF = min(refMaxRBF, origMax);
% 
% N_infill_RBF    = RBFsurrogateNodes;  % 5 new nodes inside the refined space
% Dimensions      = 3;
% numSets_inf_RBF = 2000;
% 
% bestSD_inf_RBF  = inf; %Dstar!
% bestXn_inf_RBF  = [];
% 
% for k = 1:numSets_inf_RBF
%     Xn_cand_inf_RBF = rand(N_infill_RBF, Dimensions);   % in [0,1]^3
%     sD_cand_inf_RBF = starD(Xn_cand_inf_RBF);           % star discrepancy
% 
%     if sD_cand_inf_RBF < bestSD_inf_RBF
%         bestSD_inf_RBF  = sD_cand_inf_RBF;
%         bestXn_inf_RBF  = Xn_cand_inf_RBF;
%     end
% end
% 
% samplesInsideWindow_RBF = zeros(N_infill_RBF, Dimensions); 
% samplesInsideWindow_RBF(:,1) = refMinRBF(1) + bestXn_inf_RBF(:,1) * (refMaxRBF(1) - refMinRBF(1));
% samplesInsideWindow_RBF(:,2) = refMinRBF(2) + bestXn_inf_RBF(:,2) * (refMaxRBF(2) - refMinRBF(2));
% samplesInsideWindow_RBF(:,3) = refMinRBF(3) + bestXn_inf_RBF(:,3) * (refMaxRBF(3) - refMinRBF(3));
% 
% MCSI_infill_RBF = ChosenFunc(samplesInsideWindow_RBF); %Evaluate the values of new smples
% 
% samples_RBF_new      = [samples; samplesInsideWindow_RBF]; %Add old values with the new
% MCSI_samples_RBF_new = [MCSI_samples; MCSI_infill_RBF];
% 
% N_RBF_new = size(samples_RBF_new, 1);
% 
% X_norm_RBF_new = (samples_RBF_new - Xmin) ./ (Xmax - Xmin); % Normalise with same global Xmin, Xmax
% 
% D_new = norm(max(X_norm_RBF_new) - min(X_norm_RBF_new)); % diameter and length scale using augmented data
% L_new = D_new / (0.8 * (sqrt(N_RBF_new))^(1/4));
% 
% Phi_RBF_new = zeros(N_RBF_new, N_RBF_new);% Build new Phi_RBF_new (N_RBF_new × N_RBF_new)
% for i = 1:N_RBF_new
%     for j = 1:N_RBF_new
%         r_ij = sqrt( sum( (X_norm_RBF_new(i,:) - X_norm_RBF_new(j,:)).^2 , 2) );
%         Phi_RBF_new(i,j) = exp( - (r_ij / L_new)^2 );
%     end
% end
% 
% W_RBF_new = Phi_RBF_new \ MCSI_samples_RBF_new; % New weights W_RBF_new
% 
% Xq_RBF_ref = [ ...
%     rand(2000,1)*(maxAdvert-minAdvert) + minAdvert, ...
%     rand(2000,1)*(maxTarget-minTarget) + minTarget, ...
%     rand(2000,1)*(maxSkew-minSkew)     + minSkew   ];
% 
% Xq_norm_RBF_ref = (Xq_RBF_ref - Xmin) ./ (Xmax - Xmin);
% 
% Nq_ref = size(Xq_norm_RBF_ref,1);
% yhat_RBF_ref = zeros(Nq_ref,1);
% 
% for i = 1:Nq_ref
%     r_i = sqrt( sum( (X_norm_RBF_new - Xq_norm_RBF_ref(i,:)).^2 , 2) );
%     psi_i = exp( - (r_i / L_new).^2 );
%     yhat_RBF_ref(i) = psi_i' * W_RBF_new;
% end
% 
% %=========RBF Second Optimize!!!=================
% 
% f_RBF_ref = @(x) ...
%     exp( - ( sqrt( sum( (X_norm_RBF_new - ...
%               ((x - Xmin) ./ (Xmax - Xmin)) ).^2 , 2) ) / L_new ).^2 )' ...
%     * W_RBF_new;
% 
% [xRBF_opt2, fRBF_opt2] = ga(@(x) f_RBF_ref(x), 3, [], [], [], [], lowerBound, upperBound, [], optsGA);
% 
% xRBF_opt2 = max(xRBF_opt2, origMin);
% xRBF_opt2 = min(xRBF_opt2, origMax);
% 
% fprintf('Refined RBF surrogate GA min at [%.3f, %.3f, %.3f], surrogate = %.4f\n', xRBF_opt2(1), xRBF_opt2(2), xRBF_opt2(3), fRBF_opt2);
% 
% zoomFactorRBF2 = 0.3;  % 50% of original box size
% 
% halfWidthRBF2 = 0.5 * zoomFactorRBF2 .* (origMax - origMin);
% 
% refMinRBF2 = xRBF_opt2 - halfWidthRBF2;
% refMaxRBF2 = xRBF_opt2 + halfWidthRBF2;
% 
% refMinRBF2 = max(refMinRBF2, origMin);
% refMaxRBF2 = min(refMaxRBF2, origMax);
% 
% 
% N_infill_RBF2    = RBFthirdSurrogateNodes;   % 5 new nodes
% Dimensions       = 3;
% numSets_inf_RBF2 = 2000;
% 
% bestSD_inf_RBF2  = inf;
% bestXn_inf_RBF2  = [];
% 
% for k = 1:numSets_inf_RBF2
%     Xn_cand_inf_RBF2 = rand(N_infill_RBF2, Dimensions);  
%     sD_cand_inf_RBF2 = starD(Xn_cand_inf_RBF2);
%     if sD_cand_inf_RBF2 < bestSD_inf_RBF2
%         bestSD_inf_RBF2 = sD_cand_inf_RBF2;
%         bestXn_inf_RBF2 = Xn_cand_inf_RBF2;
%     end
% end
% 
% samplesInsideWindow_RBF2 = zeros(N_infill_RBF2, Dimensions);
% samplesInsideWindow_RBF2(:,1) = refMinRBF2(1) + bestXn_inf_RBF2(:,1) * (refMaxRBF2(1) - refMinRBF2(1));
% samplesInsideWindow_RBF2(:,2) = refMinRBF2(2) + bestXn_inf_RBF2(:,2) * (refMaxRBF2(2) - refMinRBF2(2));
% samplesInsideWindow_RBF2(:,3) = refMinRBF2(3) + bestXn_inf_RBF2(:,3) * (refMaxRBF2(3) - refMinRBF2(3));
% 
% MCSI_infill_RBF2 = ChosenFunc(samplesInsideWindow_RBF2);
% 
% % Build THIRD RBF surrogate on All combined data 
% samples_RBF_third      = [samples_RBF_new; samplesInsideWindow_RBF2];
% MCSI_samples_RBF_third = [MCSI_samples_RBF_new; MCSI_infill_RBF2];
% 
% N_RBF_third = size(samples_RBF_third, 1);
% 
% X_norm_RBF_third = (samples_RBF_third - Xmin) ./ (Xmax - Xmin);
% 
% D_third = norm(max(X_norm_RBF_third) - min(X_norm_RBF_third));
% L_third = D_third / (0.8 * (sqrt(N_RBF_third))^(1/4));
% 
% Phi_RBF_third = zeros(N_RBF_third, N_RBF_third);
% for i = 1:N_RBF_third
%     for j = 1:N_RBF_third
%         r_ij = sqrt( sum( (X_norm_RBF_third(i,:) - X_norm_RBF_third(j,:)).^2 , 2) );
%         Phi_RBF_third(i,j) = exp( - (r_ij / L_third)^2 );
%     end
% end
% 
% W_RBF_third = Phi_RBF_third \ MCSI_samples_RBF_third;
% 
% Xq_RBF_third = [ ...
%     rand(2000,1)*(maxAdvert-minAdvert) + minAdvert, ...
%     rand(2000,1)*(maxTarget-minTarget) + minTarget, ...
%     rand(2000,1)*(maxSkew-minSkew)     + minSkew   ];
% 
% Xq_norm_RBF_third = (Xq_RBF_third - Xmin) ./ (Xmax - Xmin);
% 
% Nq_third = size(Xq_norm_RBF_third,1);
% 
% yhat_RBF_third = zeros(Nq_third,1);
% for i = 1:Nq_third
%     r_i = sqrt( sum( (X_norm_RBF_third - Xq_norm_RBF_third(i,:)).^2 , 2) );
%     psi_i = exp( - (r_i / L_third).^2 );
%     yhat_RBF_third(i) = psi_i' * W_RBF_third;
% end
% 
% %======================Quality of Fit for RBF===================================
% 
% %Quality of Fit for the Third Surrogate
% yhatLOO_RBF = zeros(N_RBF_third,1);
% sseLOO_RBF  = 0;
% 
% for i = 1:N_RBF_third
%     idx_train = setdiff(1:N_RBF_third, i);
%     Phi_train = Phi_RBF_third(idx_train,:);
%     y_train   = MCSI_samples_RBF_third(idx_train);
% 
%     A_i = pinv(Phi_train) * y_train;
%     yhat_i = Phi_RBF_third(i,:) * A_i;
% 
%     yhatLOO_RBF(i) = yhat_i;
%     sseLOO_RBF = sseLOO_RBF + (MCSI_samples_RBF_third(i) - yhat_i)^2;
% end
% 
% RMSE_LOO_third = sqrt(sseLOO_RBF / N_RBF_third);
% minM3 = min(MCSI_samples_RBF_third);
% maxM3 = max(MCSI_samples_RBF_third);
% rangeM3 = maxM3 - minM3;
% 
% errorPercent_third = (RMSE_LOO_third / rangeM3)*100;
% 
% fprintf('LOOCV Error For THIRD RBF Refined Surrogate = %.4f (%.2f%% of MCSI range)\n',RMSE_LOO_third, errorPercent_third);
% 
% 
% %================Find the Minima based on GA and FUll surrogate=================================
% f_RBF_ref2 = @(x) ...
%     exp( - ( sqrt( sum( (X_norm_RBF_third - ...
%               ((x - Xmin) ./ (Xmax - Xmin)) ).^2 , 2) ) / L_third ).^2 )' ...
%     * W_RBF_third;
% 
% [xRBF_opt3, fRBF_opt3] = ga(@(x) f_RBF_ref2(x), 3, [], [], [], [], lowerBound, upperBound, [], optsGA);
% 
% xRBF_opt3 = max(xRBF_opt3, origMin);
% xRBF_opt3 = min(xRBF_opt3, origMax);
% 
% fprintf('Second Refined RBF surrogate GA min at [%.3f, %.3f, %.3f], surrogate = %.4f\n', xRBF_opt3(1), xRBF_opt3(2), xRBF_opt3(3), fRBF_opt3);
% 
% figure('Name','Third RBF Surrogate (Second Refinement)','NumberTitle','off');
% scatter3(Xq_RBF_third(:,1), Xq_RBF_third(:,2), Xq_RBF_third(:,3), 25, yhat_RBF_third, 'filled');
% hold on;
% 
% %====================Draw refinement window=====================
% 
% x1_minRBF = refMinRBF(1); x1_maxRBF = refMaxRBF(1);
% x2_minRBF = refMinRBF(2); x2_maxRBF = refMaxRBF(2);
% x3_minRBF = refMinRBF(3); x3_maxRBF = refMaxRBF(3);
% 
% % 8 corner points of the box
% C = [ x1_minRBF x2_minRBF x3_minRBF;  % 1
%       x1_maxRBF x2_minRBF x3_minRBF;  % 2
%       x1_maxRBF x2_maxRBF x3_minRBF;  % 3
%       x1_minRBF x2_maxRBF x3_minRBF;  % 4
%       x1_minRBF x2_minRBF x3_maxRBF;  % 5
%       x1_maxRBF x2_minRBF x3_maxRBF;  % 6
%       x1_maxRBF x2_maxRBF x3_maxRBF;  % 7
%       x1_minRBF x2_maxRBF x3_maxRBF]; % 8
% 
% % Faces of the box (each row indexes into C)
% faces = [1 2 3 4;  % bottom (z = min)
%          5 6 7 8;  % top    (z = max)
%          1 2 6 5;  % side x+
%          2 3 7 6;  % side y+
%          3 4 8 7;  % side x-
%          4 1 5 8]; % side y-
% 
% % Plot semi-transparent refinement window
% patch('Vertices', C, 'Faces', faces, ...
%       'FaceColor', 'none', ...        % no fill, just edges OR:
%       'EdgeColor', 'r', ...
%       'LineWidth', 1.5);
% 
% x1_min3RBF = refMinRBF2(1); x1_max3RBF = refMaxRBF2(1);
% x2_min3RBF = refMinRBF2(2); x2_max3RBF = refMaxRBF2(2);
% x3_min3RBF = refMinRBF2(3); x3_max3RBF = refMaxRBF2(3);
% 
% % 8 corner points of the box
% C = [ x1_min3RBF x2_min3RBF x3_min3RBF;  % 1
%       x1_max3RBF x2_min3RBF x3_min3RBF;  % 2
%       x1_max3RBF x2_max3RBF x3_min3RBF;  % 3
%       x1_min3RBF x2_max3RBF x3_min3RBF;  % 4
%       x1_min3RBF x2_min3RBF x3_max3RBF;  % 5
%       x1_max3RBF x2_min3RBF x3_max3RBF;  % 6
%       x1_max3RBF x2_max3RBF x3_max3RBF;  % 7
%       x1_min3RBF x2_max3RBF x3_max3RBF]; % 8
% 
% % Faces of the box (each row indexes into C)
% faces3 = [1 2 3 4;  % bottom (z = min)
%          5 6 7 8;  % top    (z = max)
%          1 2 6 5;  % side x+
%          2 3 7 6;  % side y+
%          3 4 8 7;  % side x-
%          4 1 5 8]; % side y-
% 
% % Plot semi-transparent refinement window
% patch('Vertices', C, 'Faces', faces3, ...
%       'FaceColor', 'none', ...        % no fill, just edges OR:
%       'EdgeColor', 'g', ...
%       'LineWidth', 1.5);
% 
% %==========================================================
% 
% xlabel('x1 (Advert)');
% ylabel('x2 (Target)');
% zlabel('x3 (Skew)');
% title('RBF 4D Scatter of Third RBF Surrogate (colour = predicted MCSI)');
% colorbar;
% grid on;
% view(45,25);
% hold off;
% 
% %===============Plot the Original MCSI samples VS the NEW=============================
% figure('Name', 'RBF Original Samples Vs Refined1 and Refined2', 'NumberTitle', 'off');
% hold on;
% 
% % --- Original samples (all initial star discrepancy points) ---
% scatter3(samples(:,1), samples(:,2), samples(:,3), ...
%          60, MCSI_samples, 'o', ...
%          'MarkerFaceColor', [0 0.447 0.741], ...   % blue
%          'MarkerEdgeColor', 'k', ...
%          'DisplayName', 'Original Samples');
% 
% % --- New infill samples (inside window) ---
% scatter3(samplesInsideWindow_RBF(:,1), samplesInsideWindow_RBF(:,2), samplesInsideWindow_RBF(:,3), ...
%          80, MCSI_infill, 's', ...
%          'MarkerFaceColor', [0.85 0.325 0.098], ... % orange/red
%          'MarkerEdgeColor', 'k', ...
%          'DisplayName', 'Infill Samples');
% 
% % --- New infill samples (inside window) ---
% scatter3(samplesInsideWindow_RBF2(:,1), samplesInsideWindow_RBF2(:,2), samplesInsideWindow_RBF2(:,3), ...
%          80, MCSI_infill, 's', ...
%          'MarkerFaceColor', [0.4660 0.6740 0.1880], ... % orange/red
%          'MarkerEdgeColor', 'k', ...
%          'DisplayName', 'Infill Samples');
% 
% xlabel('x1 (Advert)');
% ylabel('x2 (Target)');
% zlabel('x3 (Skew)');
% title('RBF Original vs Infill Sample Points (coloured by MCSI)');
% colorbar;
% grid on;
% view(45,25);
% legend('Location','best');
% hold off;
% 
% %FInd the True Minima for the example functions
% [xTrue, fTrue] = ga(@(x) ChosenFunc(x), 3, [], [], [], [], lowerBound, upperBound, [], optsGA);
% 
% fprintf('GA Function True min at [%.3f, %.3f, %.3f]\n', xTrue(1), xTrue(2), xTrue(3));
% fprintf('  True MCSI  = %.4f\n', fTrue);
