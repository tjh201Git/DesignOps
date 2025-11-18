% draw refinement box
function drawRefinementBox(refMin,refMax, colour)
% Draw refinement box
        % Vertices of the cube
        xBox = [refMin(1) refMax(1) refMax(1) refMin(1) refMin(1) refMax(1) refMax(1) refMin(1)];
        yBox = [refMin(2) refMin(2) refMax(2) refMax(2) refMin(2) refMin(2) refMax(2) refMax(2)];
        zBox = [refMin(3) refMin(3) refMin(3) refMin(3) refMax(3) refMax(3) refMax(3) refMax(3)];
        
        % Define faces
        faces = [1 2 3 4;
                 5 6 7 8;
                 1 2 6 5;
                 2 3 7 6;
                 3 4 8 7;
                 4 1 5 8];
        
        patch('Vertices',[xBox', yBox', zBox'], 'Faces',faces, ...
              'FaceColor','none','EdgeColor',colour,'LineWidth',2,'LineStyle','--');

end