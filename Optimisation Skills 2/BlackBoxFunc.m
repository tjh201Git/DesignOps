
advertSpend = 1%Input 1, 0-200
targetSpend = 2%Input 2, 0-60
audienceSkew = 3 %Input 3, 0-1

%first select 3 values as inputs
%input 3 values into a blackbox function
%plot the MCSI result
%based on the method and the MCSI the next 3 a chosen
%repeat until the highest MCSI is found

%Matlab GA function



fsurf(@(x, y) reshape(exampFunc([x(:), y(:)]), size(x)), ...
    'MeshDensity', 100, ...
    'ShowContours', 'on', ...
    'LineStyle', ':');

xlabel('x');
ylabel('y');
zlabel('exampFunc(x, y)');
title('Surface Plot of exampFunc');