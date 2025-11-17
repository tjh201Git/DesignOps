% find the best distribution of points with the lowest D* val

function [pointsMatrix, bestDStar] = bestDStarPointDistribution(dimensions, nodes, iterations, saveFile)

    if (nargin < 4)
        saveFile = false;
    end
    
    bestDStar  = inf; %Start with infinity then replace with smaller and smaller discrepency values
    pointsMatrix  = []; %placeholder for best design
    
    for k = 1:iterations %For each set
    
        Xn_cand = rand(nodes, dimensions); %Generate a random N×D point set in the hypercube
        sD_cand = starD(Xn_cand); %Find the star disrepency
    
        if sD_cand < bestDStar %Keep the smallest discrepency so far
            bestDStar = sD_cand;
            pointsMatrix = Xn_cand; 
        end
    end
    
    fprintf('Best star discrepancy (normalised space) = %.4f\n', bestDStar);

    if saveFile
        distribution.bestDStar = bestDStar;
        distribution.pointsMatrix = pointsMatrix;

        fileName = "points_" + dimensions + "dims_" + nodes + "nodes.mat";

        save(fileName,"distribution");


    end

end