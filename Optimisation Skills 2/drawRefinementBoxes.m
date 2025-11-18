function drawRefinementBoxes(refMins, refMaxs)
% DRAwREFINEMENTBOXES Draws a series of refinement boxes using the parula colormap.
% refMins: 2xN matrix of min coordinates [x_min; y_min] for N boxes.
% refMaxs: 2xN matrix of max coordinates [x_max; y_max] for N boxes.
    
    % 1. Determine the number of refinement boxes
    numRefs = size(refMins, 2); 
    
    if numRefs == 0
        return; % Nothing to draw
    end
    
    % 2. Generate the parula color array (N x 3)
    % colormap(name) generates the colormap, and colormap(name, N) sets the size.
    % Using parula(N) directly returns an N x 3 array of RGB values.
    colors = winter(numRefs); 
    
    % 3. Loop through and draw each box
    for i = 1:numRefs
        % Get the color for the current box (i-th row of the colors array)
        currentColor = colors(i, :); 

        refMin = refMins(i, :);
        refMax = refMaxs(i, :);
        
        % Pass the color and coordinates to the helper function
        drawRefinementBox(refMin, refMax, currentColor);
    end
    
end