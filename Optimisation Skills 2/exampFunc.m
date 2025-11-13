function output = exampFunc(x)
    output = 10.0 * size(x,2) + sum(x .^2 - 10.0 * cos(2 * pi .* x),2);
