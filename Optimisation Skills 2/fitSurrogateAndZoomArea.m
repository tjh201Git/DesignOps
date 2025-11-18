function [refMin, refMax] = fitSurrogateAndZoomArea(samples, MCSI_samples, zoomFactor, plotGraphs, stageName)

    if nargin < 4
        stageName = string(length(samples(:,1)));
    end
    if nargin < 3
        plotGraphs = true;
    end

    % ================================================
    % ====== Define bounds and normalize ============
    % ================================================
    minAdvert = 0; maxAdvert = 200;
    minTarget = 0; maxTarget = 60;
    minSkew = 0; maxSkew = 1;

    Xmin = [minAdvert, minTarget, minSkew];
    Xmax = [maxAdvert, maxTarget, maxSkew];

    X_norm = (samples - Xmin) ./ (Xmax - Xmin);
    Nodes = length(samples(:,1));

    x1 = X_norm(:,1); x2 = X_norm(:,2); x3 = X_norm(:,3);

    optsGA = optimoptions('ga', 'MaxGenerations', 100, 'PopulationSize', 1000);

    % ================================================
    % ======= Polynomial Surrogate ==================
    % ================================================
    Phi_orders = cell(5,1);
    RMSE_poly = zeros(5,1);
    for order = 1:5
        Phi_orders{order} = buildPhi(Nodes, order, x1, x2, x3);
        RMSE_poly(order) = LOOV(Nodes, Phi_orders{order}, MCSI_samples);
    end
    [bestRMSE_poly, BestOrder] = min(RMSE_poly);
    errorPercent = (bestRMSE_poly / (max(MCSI_samples) - min(MCSI_samples)))*100;
    Phi_poly = Phi_orders{BestOrder};
    A_poly = pinv(Phi_poly) * MCSI_samples;

    % ================================================
    % ======= RBF Surrogate =========================
    % ================================================
    epsilon = 1; % Gaussian RBF
    Phi_rbf = zeros(Nodes, Nodes);
    for i = 1:Nodes
        for j = 1:Nodes
            r = norm(X_norm(i,:) - X_norm(j,:));
            Phi_rbf(i,j) = exp(-(epsilon*r)^2);
        end
    end
    w_rbf = Phi_rbf \ MCSI_samples;

    % Vectorized LOOCV for RBF
    Phi_inv = inv(Phi_rbf);
    yhat_full = Phi_rbf * w_rbf;
    loo_residuals = (MCSI_samples - yhat_full) ./ diag(Phi_inv);
    RMSE_rbf = sqrt(mean(loo_residuals.^2));






    % % ================================================
    % % ======= Choose Best Surrogate =================
    % % ================================================
    % if bestRMSE_poly <= RMSE_rbf
    %     surrogateType = "poly";
    %     RMSE_best = bestRMSE_poly;
    %     A_best = A_poly;
    %     w_best = [];          % not used
    %     BestOrder_poly = BestOrder;
    %     X_norm_rbf = [];      % not used
    %     fprintf("Using Polynomial Surrogate (Order %d), RMSE_LOO=%.4f\n", BestOrder, bestRMSE_poly);
    % else
    %     surrogateType = "rbf";
    %     RMSE_best = RMSE_rbf;
    %     w_best = w_rbf;
    %     A_best = [];          % not used
    %     BestOrder_poly = [];  % not used
    %     X_norm_rbf = X_norm;  % needed for GA evaluation
    %     fprintf("Using RBF Surrogate, RMSE_LOO=%.4f\n", RMSE_rbf);
    % end
    % 
    % % ================================================
    % % ======= Build Query Points ====================
    % % ================================================
    % Nq = 2000;
    % Xq = [ ...
    %     rand(Nq,1)*(maxAdvert-minAdvert)+minAdvert, ...
    %     rand(Nq,1)*(maxTarget-minTarget)+minTarget, ...
    %     rand(Nq,1)*(maxSkew-minSkew)+minSkew];
    % Xq_norm = (Xq - Xmin) ./ (Xmax - Xmin);
    % x1qn = Xq_norm(:,1); x2qn = Xq_norm(:,2); x3qn = Xq_norm(:,3);
    % 
    % % Evaluate surrogate at query points
    % if surrogateType == "poly"
    %     Phi_q = buildPhi_q(Xq, BestOrder_poly, x1qn, x2qn, x3qn);
    %     yhat = Phi_q * A_best;
    % else
    %     % RBF evaluation
    %     Phi_q = zeros(Nq, Nodes);
    %     for i = 1:Nq
    %         for j = 1:Nodes
    %             r = norm(Xq_norm(i,:) - X_norm_rbf(j,:));
    %             Phi_q(i,j) = exp(-(epsilon*r)^2);
    %         end
    %     end
    %     yhat = Phi_q * w_best;
    % end




    % ================================================
    % ======= Build Query Points ====================
    % ================================================
    Nq = 2000;
    Xq = [ ...
        rand(Nq,1)*(maxAdvert-minAdvert)+minAdvert, ...
        rand(Nq,1)*(maxTarget-minTarget)+minTarget, ...
        rand(Nq,1)*(maxSkew-minSkew)+minSkew];
    Xq_norm = (Xq - Xmin) ./ (Xmax - Xmin);
    x1qn = Xq_norm(:,1); x2qn = Xq_norm(:,2); x3qn = Xq_norm(:,3);

    % Evaluate surrogate at query points: poly
        Phi_q_poly = buildPhi_q(Xq, BestOrder, x1qn, x2qn, x3qn);
        yhat_poly = Phi_q_poly * A_poly;

        % RBF evaluation
        Phi_q_rbf = zeros(Nq, Nodes);
        for i = 1:Nq
            for j = 1:Nodes
                r = norm(Xq_norm(i,:) - X_norm(j,:));
                Phi_q_rbf(i,j) = exp(-(epsilon*r)^2);
            end
        end
        yhat_rbf = Phi_q_rbf * w_rbf;

        
    % ================================================
    % ======= GA Optimization Poly =======================
    % ================================================
    f_sur = @(x) surrogateEval(x, 'poly', Xmin, Xmax, Nodes, A_poly, w_rbf, BestOrder, X_norm, epsilon);

    lowerBound = Xmin; upperBound = Xmax;
    [xSur_opt_poly, f_opt_poly] = ga(f_sur, 3, [], [], [], [], lowerBound, upperBound, [], optsGA);

    % Clip within original bounds
    xSur_opt_poly = max(xSur_opt_poly, Xmin);
    xSur_opt_poly = min(xSur_opt_poly, Xmax);
    fprintf("Best guess according to surrogate: MSCI = %.4f @ (%.4f %.4f %.4f)\n", f_opt_poly, xSur_opt_poly);

    % ================================================
    % ======= Zoom Refinement =======================
    % ================================================
    % zoomFactor = 0.5;
    halfWidth = 0.5 * zoomFactor .* (Xmax - Xmin);
    refMin_poly = max(xSur_opt_poly - halfWidth, Xmin);
    refMax_poly = min(xSur_opt_poly + halfWidth, Xmax);

    % ================================================
    % ======= GA Optimization RBF=======================
    % ================================================
    f_sur = @(x) surrogateEval(x, 'rbf', Xmin, Xmax, Nodes, A_poly, w_rbf, BestOrder, X_norm, epsilon);

    lowerBound = Xmin; upperBound = Xmax;
    [xSur_opt_rbf, f_opt_rbf] = ga(f_sur, 3, [], [], [], [], lowerBound, upperBound, [], optsGA);

    % Clip within original bounds
    xSur_opt_rbf = max(xSur_opt_rbf, Xmin);
    xSur_opt_rbf = min(xSur_opt_rbf, Xmax);
    fprintf("Best guess according to surrogate: MSCI = %.4f @ (%.4f %.4f %.4f)\n", f_opt_rbf, xSur_opt_rbf);

    % ================================================
    % ======= Zoom Refinement =======================
    % ================================================
    % zoomFactor = 0.5;
    halfWidth = 0.5 * zoomFactor .* (Xmax - Xmin);
    refMin_rbf = max(xSur_opt_poly - halfWidth, Xmin);
    refMax_rbf = min(xSur_opt_poly + halfWidth, Xmax);
    
    switch stageName
    case "30 node Exploration"
        boxColor = 'r';
    case "Additional 5 Nodes in Refined Zone"
        boxColor = 'r';
    case "Additional 5 Nodes in Further Refined Zone"
        boxColor = 'g';
    end
    

    nexttile; %poly tile
    scatter3(Xq(:,1), Xq(:,2), Xq(:,3), 25, yhat_poly, 'filled'); 
    hold on;

    plot3(xSur_opt_poly(1), xSur_opt_poly(2), xSur_opt_poly(3), ...
    'rp', 'MarkerSize', 18, 'MarkerFaceColor', 'y', 'DisplayName','Poly Minima');


    minLabel = sprintf('Minima:(%.2f, %.2f, %.2f), Best Order:%f LOOV Percentage Error: %.1f%%,', ...
           xSur_opt_poly(1), xSur_opt_poly(2), xSur_opt_poly(3), BestOrder, errorPercent);           

    text(xSur_opt_poly(1), xSur_opt_poly(2), xSur_opt_poly(3), ...
    minLabel, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom', ...
    'FontWeight','bold');
    hold on;

    xlabel('x1 (Advert)'); 
    ylabel('x2 (Target)'); 
    zlabel('x3 (Skew)');
    title(sprintf('%s Surrogate Model, Method: %s', stageName, 'poly'));
    colorbar; 
    grid on; 
    view(45,25);

    xlim([0 200]);
    ylim([0 60])
    zlim([0 1]);

    drawRefinementBox(refMin_poly, refMax_poly, boxColor);

    nexttile; %rbf tile
    scatter3(Xq(:,1), Xq(:,2), Xq(:,3), 25, yhat_rbf, 'filled'); 
    hold on;

    plot3(xSur_opt_rbf(1), xSur_opt_rbf(2), xSur_opt_rbf(3), ...
    'rp', 'MarkerSize', 18, 'MarkerFaceColor', 'y', 'DisplayName','Poly Minima');


    minLabel = sprintf('Minima:(%.2f, %.2f, %.2f)', ...
           xSur_opt_rbf(1), xSur_opt_rbf(2), xSur_opt_rbf(3));           

    text(xSur_opt_rbf(1), xSur_opt_rbf(2), xSur_opt_rbf(3), ...
    minLabel, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom', ...
    'FontWeight','bold');
    hold on;

    xlabel('x1 (Advert)'); 
    ylabel('x2 (Target)'); 
    zlabel('x3 (Skew)');
    title(sprintf('%s Surrogate Model, Method: %s', stageName, 'rbf'));
    colorbar; 
    grid on; 
    view(45,25);

    xlim([0 200]);
    ylim([0 60])
    zlim([0 1]);

    drawRefinementBox(refMin_rbf, refMax_rbf, boxColor);

    

     while true %wait for the user to input MCSI values
        userInput = input('Enter if you like "poly" or "rbf" more: ', 's');
        
        if strcmpi(userInput, 'poly')
            surrogateType = "poly";
            break;   % Exit loop
        elseif strcmpi(userInput, 'rbf')
            surrogateType = "rbf";
            break;   % Exit loop
        else
            pause(1)

        end 

     end
% 
% surrogateType = "poly";
%         RMSE_best = bestRMSE_poly;
%         A_best = A_poly;
%         w_best = [];          % not used
%         BestOrder_poly = BestOrder;
%         X_norm_rbf = [];      % not used
%         fprintf("Using Polynomial Surrogate (Order %d), RMSE_LOO=%.4f\n", BestOrder, bestRMSE_poly);
   % ================================================
    % ======= Choose Best Surrogate =================
    % ================================================
    if userInput == "poly"
        % plotSurrogateAndSampleGraph(Xq, yhat_poly, samples, 'poly', xSur_opt_poly, refMin_poly, refMax_poly, surrogateType, errorPercent, BestOrder);
    else
        errorPercent = 0;
        BestOrder = 0;
        % plotSurrogateAndSampleGraph(Xq, yhat_rbf, samples, 'rbf', xSur_opt_rbf, refMin_rbf, refMax_rbf, surrogateType, errorPercent, BestOrder);

    end


    % ================================================
    % ======= Optional Plot =========================
    % ================================================
    % if plotGraphs
    %     figure('Name', 'Surrogate Model', 'NumberTitle', 'off');
    % 
    %     % Surrogate predicted points
    %     scatter3(Xq(:,1), Xq(:,2), Xq(:,3), 25, yhat, 'filled'); 
    %     hold on;
    % 
    %     % Actual sample points in white
    %     scatter3(samples(:,1), samples(:,2), samples(:,3), 60, 'w', 'filled', 'MarkerEdgeColor','k');
    % 
    %     xlabel('x1 (Advert)'); 
    %     ylabel('x2 (Target)'); 
    %     zlabel('x3 (Skew)');
    %     title('Surrogate Model Predictions (cube) with Actual Samples');
    %     colorbar; 
    %     grid on; 
    %     view(45,25);
    % end

    % ================================================
    % ======= GA Optimization =======================
    % ================================================
    % f_sur = @(x) surrogateEval(x, surrogateType, Xmin, Xmax, Nodes, A_poly, w_rbf, BestOrder, X_norm, epsilon);
    % 
    % lowerBound = Xmin; upperBound = Xmax;
    % [xSur_opt, f_opt] = ga(f_sur, 3, [], [], [], [], lowerBound, upperBound, [], optsGA);
    % 
    % % Clip within original bounds
    % xSur_opt = max(xSur_opt, Xmin);
    % xSur_opt = min(xSur_opt, Xmax);
    % fprintf("Best guess according to surrogate: MSCI = %.4f @ (%.4f %.4f %.4f)\n", f_opt, xSur_opt);
    % 
    % % ================================================
    % % ======= Zoom Refinement =======================
    % % ================================================
    % % zoomFactor = 0.5;
    % halfWidth = 0.5 * zoomFactor .* (Xmax - Xmin);
    % refMin = max(xSur_opt - halfWidth, Xmin);
    % refMax = min(xSur_opt + halfWidth, Xmax);

    % ================================================
    % ======= Optional Plot =========================
    % ================================================
    if strcmpi(surrogateType, 'poly')
        yhat = yhat_poly;
        refMin = refMin_poly;
        refMax = refMax_poly;
        xSur_opt = xSur_opt_poly;
    else
        yhat = yhat_rbf;
        refMin = refMin_rbf;
        refMax = refMax_rbf;
        xSur_opt = xSur_opt_rbf;
    end

    if plotGraphs
        plotSurrogateAndSampleGraph(Xq, yhat, samples, stageName, xSur_opt, refMin, refMax, surrogateType, errorPercent, BestOrder);
    end


end

%% ===== Helper Function for GA Evaluation =====
function y = surrogateEval(x, type, Xmin, Xmax, Nodes, A_poly, w_rbf, BestOrder, X_norm_rbf, epsilon, Phi_order)
    if nargin < 10
        epsilon = 1;
    end
    x_norm = (x - Xmin) ./ (Xmax - Xmin);
    if type == "poly"
        Phi_x = buildPhi_q(x_norm, BestOrder, x_norm(1), x_norm(2), x_norm(3));
        % Phi_x = [1, x_norm(1), x_norm(2), x_norm(3), ...
        %          x_norm(1)^2, x_norm(2)^2, x_norm(3)^2, ...
        %          x_norm(1)*x_norm(2), x_norm(1)*x_norm(3), x_norm(2)*x_norm(3)];
        y = Phi_x * A_poly;
    else
        % RBF evaluation
        Phi_x = zeros(1, Nodes);
        for j = 1:Nodes
            r = norm(x_norm - X_norm_rbf(j,:));
            Phi_x(j) = exp(-(epsilon*r)^2);
        end
        y = Phi_x * w_rbf;
    end
end
