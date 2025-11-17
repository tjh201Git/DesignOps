function output = ackley3(x)
    x1 = x(:,1); x2 = x(:,2); x3 = x(:,3);

    term1 = -20 * exp(-0.2 * sqrt((x1.^2 + x2.^2 + x3.^2) / 3));
    term2 = -exp((cos(2*pi*x1) + cos(2*pi*x2) + cos(2*pi*x3)) / 3);

    output = term1 + term2 + 20 + exp(1);
end
