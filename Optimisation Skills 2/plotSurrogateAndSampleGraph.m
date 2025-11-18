% plot surrogate and sample points
function plotSurrogateAndSampleGraph(surrogatePoints, yhat, samplePoints)
        figure('Name', 'Surrogate Model', 'NumberTitle', 'off');
        
        % Surrogate predicted points (cube)
        scatter3(surrogatePoints(:,1), surrogatePoints(:,2), surrogatePoints(:,3), 25, yhat, 'filled'); 
        hold on;
        
        % Actual sample points in white
        scatter3(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3), 60, 'w', 'filled', 'MarkerEdgeColor','k');
        
        
        xlabel('x1 (Advert)'); 
        ylabel('x2 (Target)'); 
        zlabel('x3 (Skew)');
        title('Surrogate Model with Actual Samples and Refinement Box');
        colorbar; 
        grid on; 
        view(45,25);

end