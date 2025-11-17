function output = griewank3(x)
    x1 = x(:,1); x2 = x(:,2); x3 = x(:,3);

    sumTerm = (x1.^2 + x2.^2 + x3.^2) / 4000;
    prodTerm = cos(x1./1) .* cos(x2./sqrt(2)) .* cos(x3./sqrt(3));

    output = 1 + sumTerm - prodTerm;
end
