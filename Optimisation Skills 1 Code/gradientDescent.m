function [x, fval, history] = gradientDescent(f, gradf, x0, alpha, tol, maxIter)
% GRADIENTDESCENT Generic gradient descent optimizer
%
% Inputs:
%   f       - function handle, f(x)
%   gradf   - gradient function handle, gradf(x)
%   x0      - initial guess (scalar or vector)
%   alpha   - learning rate (step size)
%   tol     - tolerance on gradient norm
%   maxIter - maximum iterations
%
% Outputs:
%   x       - final value of x
%   fval    - f(x) at solution
%   history - struct containing iteration information

    if nargin < 6, maxIter = 1000; end
    if nargin < 5, tol = 1e-6; end

    x = x0;
    history.x = zeros(length(x0), maxIter);
    history.fval = zeros(1, maxIter);
    history.gradNorm = zeros(1, maxIter);

    for k = 1:maxIter
        g = gradf(x);
        history.x(:, k) = x;
        history.fval(k) = f(x);
        history.gradNorm(k) = norm(g);

        % stopping criterion
        if norm(g) < tol
            break;
        end

        % update rule
        x = x - alpha * g;
    end

    % truncate history
    history.x = history.x(:, 1:k);
    history.fval = history.fval(1:k);
    history.gradNorm = history.gradNorm(1:k);

    fval = f(x);
end