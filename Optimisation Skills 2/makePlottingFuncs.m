% % all our plotting functions in one place. Hardcoding is chill here.
% 
% function plottingFuncs = makePlottingFuncs()
%     plottingFuncs.blackBoxPlotting = @blackBoxPlotting;
% end
% 
% 
% 
% 
% function blackBoxPlotting(resolution)
%     %Now plotting the example blackbox function with the 3 variable ranges so
%     %that we can compare agaisnt the surrogate model 
%     %its being plotted as a point cloud because its a 4 dimensional plot
% 
%     Nblackbox = 2000; %Resolution for blackbox plotting
% 
%     X_Blackbox = [ ... %%
%         denormaliseMatrix(rand(Nblackbox,1)*(maxAdvert-minAdvert) + minAdvert, ...  % x1: advert
%         rand(Nblackbox,1)*(maxTarget-minTarget) + minTarget, ...  % x2: target
%         rand(Nblackbox,1)*(maxSkew-minSkew)     + minSkew   ];    % x3: skew
% 
%     MCSI_true = ChosenFunc(X_Blackbox);
% 
%     figure('Name', 'Plot of the Example BlackBox', 'NumberTitle', 'off');
%     scatter3(X_Blackbox(:,1), X_Blackbox(:,2), X_Blackbox(:,3), 25, MCSI_true, 'filled');
%     hold on;
% 
% 
%     xlabel('Advert Spending');
%     ylabel('Target Spending');
%     zlabel('Audience Skew');
%     title('4D Plot of the example function');
%     colorbar;
%     grid on;
%     view(45, 25);
% 
%     %scatter3(samples(:,1),samples(:,2), MCSI_samples, 50, MCSI_samples,'red', 'filled');
% 
%     hold off;
% end