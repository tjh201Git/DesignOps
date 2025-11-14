%FROM LECTURE NOTES

function V = volHyperCube(p)
%calculates volume of hypercube anchored at origin
%p is corner point

 V = prod(p);
%area of hypercube is the product of coordinates of p
%simple to calculate, as hypercube is anchored at origin (0,0,0)
end