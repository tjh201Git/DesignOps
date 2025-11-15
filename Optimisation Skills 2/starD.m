%FROM LECTURE NOTES


function sD = starD(X)
%calculation of star discrepancy in N dimensions
%X should be supplied as normalised coordinates (0 to 1)
 d = zeros([size(X,1), 1]); %initialise local discrepancy vector

 %algorithm defines hypercubes with corner point, p, defined by the
 %sample points. This greatly improves the speed of this algorithm
 %the alternative is to try every possible corner point, p, in the
 %design space, which can be expensive

 for i = 1:size(X,1) %iterate through all sample points
    p = X(i,:);
%define corners of hypercube as sample point i - NB: hypercube is anchored at origin
    in = inHyperCube(X,p); %find no of points in hypercube
%calculate the star discrepancy - the fraction of points inside the polygon minus the volume of the polygon
 %calc local discrepancy
    d(i) = abs( in /size(X,1) - volHyperCube(p)/1);
%as coordinates are normaised, total design space volume is 1
 end
 sD = max(d); %the star discrepancy is the max of all local star discrepancies
end

function in = inHyperCube(X,p)
 %p is corner point of hypercube

    D = size(X,2); %dimension is column size of X

 D_in = sum(X <= p, 2); 
%determine how many dimensions are within hypercube at each point

 in = sum(D_in == D);
%for a point to be completely within hyercube, all N dimensions must lie within

end
