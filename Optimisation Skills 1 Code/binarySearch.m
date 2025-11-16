function [x_mid, f_mid, iterations] = binarySearch(func, low, high, tol, maxIter)
% BINARYSEARCH  Generic binary search / bisection method
%
%   func      - function handle, e.g. @(x) x^2 - 2
%   low, high - initial interval bounds
%   tol       - tolerance on interval size
%   maxIter   - maximum number of iterations
%
% Returns:
%   x_mid     - best estimate of solution
%   f_mid     - func(x_mid)
%   iterations - number of iterations performed

    iterations = 0;

    f_low = func(low);
    f_high = func(high);

    % Ensure the function changes sign
    if f_low * f_high > 0
        error('binarySearch:sign','Function must change sign on [low, high].');
    end

    while (high - low) > tol && iterations < maxIter
        iterations = iterations + 1;

        x_mid = 0.5*(low + high);
        f_mid = func(x_mid);

        % Decide which half interval to keep
        if f_low * f_mid <= 0
            % root is in [low, mid]
            high = x_mid;
            f_high = f_mid;
        else
            % root is in [mid, high]
            low = x_mid;
            f_low = f_mid;
        end
    end

    % Final midpoint
    x_mid = 0.5*(low + high);
    f_mid = func(x_mid);
end
