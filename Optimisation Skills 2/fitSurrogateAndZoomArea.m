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

    % ================================================
    % ======= Choose Best Surrogate =================
    % ================================================
    if bestRMSE_poly <= RMSE_rbf
        surrogateType = "poly";
        RMSE_best = bestRMSE_poly;
        A_best = A_poly;
        w_best = [];          % not used
        BestOrder_poly = BestOrder;
        X_norm_rbf = [];      % not used
        fprintf("Using Polynomial Surrogate (Order %d), RMSE_LOO=%.4f\n", BestOrder, bestRMSE_poly);
    else
        surrogateType = "rbf";
        RMSE_best = RMSE_rbf;
        w_best = w_rbf;
        A_best = [];          % not used
        BestOrder_poly = [];  % not used
        X_norm_rbf = X_norm;  % needed for GA evaluation
        fprintf("Using RBF Surrogate, RMSE_LOO=%.4f\n", RMSE_rbf);
    end

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

    % Evaluate surrogate at query points
    if surrogateType == "poly"
        Phi_q = buildPhi_q(Xq, BestOrder_poly, x1qn, x2qn, x3qn);
        yhat = Phi_q * A_best;
    else
        % RBF evaluation
        Phi_q = zeros(Nq, Nodes);
        for i = 1:Nq
            for j = 1:Nodes
                r = norm(Xq_norm(i,:) - X_norm_rbf(j,:));
                Phi_q(i,j) = exp(-(epsilon*r)^2);
            end
        end
        yhat = Phi_q * w_best;
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
    f_sur = @(x) surrogateEval(x, surrogateType, Xmin, Xmax, Nodes, A_best, w_best, BestOrder_poly, X_norm_rbf, epsilon);

    lowerBound = Xmin; upperBound = Xmax;
    [xSur_opt, f_opt] = ga(f_sur, 3, [], [], [], [], lowerBound, upperBound, [], optsGA);

    % Clip within original bounds
    xSur_opt = max(xSur_opt, Xmin);
    xSur_opt = min(xSur_opt, Xmax);
    fprintf("Best guess according to surrogate: MSCI = %.4f @ (%.4f %.4f %.4f)\n", f_opt, xSur_opt);

    % ================================================
    % ======= Zoom Refinement =======================
    % ================================================
    % zoomFactor = 0.5;
    halfWidth = 0.5 * zoomFactor .* (Xmax - Xmin);
    refMin = max(xSur_opt - halfWidth, Xmin);
    refMax = min(xSur_opt + halfWidth, Xmax);

    % ================================================
    % ======= Optional Plot =========================
    % ================================================
    if plotGraphs
        plotSurrogateAndSampleGraph(Xq, yhat, samples);
        drawRefinementBox(refMin, refMax);
    end


end

%% ===== Helper Function for GA Evaluation =====
function y = surrogateEval(x, type, Xmin, Xmax, Nodes, A_poly, w_rbf, BestOrder, X_norm_rbf, epsilon)
    if nargin < 10
        epsilon = 1;
    end
    x_norm = (x - Xmin) ./ (Xmax - Xmin);
    if type == "poly"
        Phi_x = [1, x_norm(1), x_norm(2), x_norm(3), ...
                 x_norm(1)^2, x_norm(2)^2, x_norm(3)^2, ...
                 x_norm(1)*x_norm(2), x_norm(1)*x_norm(3), x_norm(2)*x_norm(3)];
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
