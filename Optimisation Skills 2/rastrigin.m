function output = rastrigin(x)
A = 10;
n = size(x,2);                   
output = A * n + sum(x.^2 - A .* cos(2*pi.*x), 2);
end
