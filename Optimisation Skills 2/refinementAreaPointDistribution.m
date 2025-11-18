function [pointsMatrix, bestDStar, numExistingSamplesInRefArea] = refinementAreaPointDistribution(samples, refMin, refMax, varMins, varMaxs, dimensions, nodes, iterations)
    
    % The 'samples' input is UNNORMALIZED data from the entire space.
    % 'refMin' and 'refMax' are the UNNORMALIZED bounds of the refinement area.
    % 'varMins' and 'varMaxs' are the UNNORMALIZED bounds of the ENTIRE problem space.

    % Initialize normWithinBoxSamples to an empty matrix for the 'no existing samples' case.
    normWithinBoxSamples = []; 

    % 1. Extract and Normalize existing samples within the UNNORMALIZED ref box
    withinBoxSamples = [];
    len = length(samples(:,1));
    for i = 1:len
        % Check is done in UNNORMALIZED space
        if isSampleWithinBounds(samples(i,:), refMin, refMax)
            withinBoxSamples = [withinBoxSamples; samples(i, :)];
        end
    end
    
    % Check if any samples were found within the refinement area
    numExistingSamplesInRefArea = size(withinBoxSamples, 1); 

    % Convert the extracted UNNORMALIZED points to the NORMALIZED space [0, 1]^D
    % This normalization uses the bounds of the entire problem space.
    if numExistingSamplesInRefArea > 0
        normWithinBoxSamples = normaliseMatrix(withinBoxSamples, varMins, varMaxs);
    end
    
    bestDStar = Inf;
    
    for i = 1:iterations
        
        % 2. Generate new random points in the refinement area
        randPoints = rand(nodes, dimensions); % Points in [0, 1]^D
        
        % Map [0, 1]^D directly into the UNNORMALIZED refinement box [refMin, refMax]
        scaledToRefBox = denormaliseMatrix(randPoints, refMin, refMax); 
        
        % 3. Convert the new UNNORMALIZED points back to the NORMALIZED space [0, 1]^D
        % This uses the bounds of the entire problem space.
        normNewPoints = normaliseMatrix(scaledToRefBox, varMins, varMaxs);
        
        % 4. Combine all points (both existing and new) in the NORMALIZED space
        % This works correctly even if normWithinBoxSamples is empty.
        normTempPointsMatrix = [normWithinBoxSamples; normNewPoints];
        
        % Optional: Sanity check to ensure points are within the [0, 1] range 
        normTempPointsMatrix(normTempPointsMatrix > 1) = 1; 
        normTempPointsMatrix(normTempPointsMatrix < 0) = 0; 
        
        % Ensure there are points to calculate D* for (only an issue if nodes=0)
        if size(normTempPointsMatrix, 1) > 0 
            dstar = starD(normTempPointsMatrix);
            if (dstar < bestDStar)
                % pointsMatrix is returned in the NORMALIZED space
                pointsMatrix = normTempPointsMatrix;
                bestDStar = dstar;
            end
        end
    
    end
    pointsMatrix = normNewPoints;

end

% --- Supporting Functions (Provided in Original Prompt) ---

function withinBounds = isSampleWithinBounds(sample, refMin, refMax)
    % Check if ALL elements of 'sample' are >= 'refMin' AND <= 'refMax'.
    % The 'all' function ensures every dimension passes the test.
    withinBounds = all(sample >= refMin) && all(sample <= refMax);
end