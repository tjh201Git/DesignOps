% plot surrogate and sample points
function plotSurrogateAndSampleGraph(surrogatePoints, yhat, samplePoints, stageName, Minima,refMin, refMax, surrogateType, errorPercent, BestOrder_poly)
        
        if (strcmpi(surrogateType, 'rbf'))
            nexttile
        end

        switch stageName
            case "30 node Exploration"
                boxColor = 'r';
            case "Additional 5 Nodes in Refined Zone"
                boxColor = 'r';
            case "Additional 5 Nodes in Further Refined Zone"
                boxColor = 'g';
        end

        %figure('Name', 'Surrogate Model', 'NumberTitle', 'off');
        % nexttile;
        % % Surrogate predicted points (cube)
        % scatter3(surrogatePoints(:,1), surrogatePoints(:,2), surrogatePoints(:,3), 25, yhat, 'filled'); 
        % hold on;
        % 
        % 
        % plot3(Minima(1), Minima(2), Minima(3), ...
        % 'rp', 'MarkerSize', 18, 'MarkerFaceColor', 'y', 'DisplayName','Poly Minima');
        % 
        % switch surrogateType
        %     case "poly"
        %     minLabel = sprintf('Minima:(%.2f, %.2f, %.2f), Best Order:%f LOOV Percentage Error: %.1f%%,', ...
        %            Minima(1), Minima(2), Minima(3), BestOrder_poly, errorPercent);           
        %     case "rbf"
        %     minLabel = sprintf('Minima:(%.2f, %.2f, %.2f)', ...
        %            Minima(1), Minima(2), Minima(3));        
        % end
        % 
        % 
        % text(Minima(1), Minima(2), Minima(3), ...
        % minLabel, ...
        % 'HorizontalAlignment','left', ...
        % 'VerticalAlignment','bottom', ...
        % 'FontWeight','bold');
        % hold on;
        % 
        % xlabel('x1 (Advert)'); 
        % ylabel('x2 (Target)'); 
        % zlabel('x3 (Skew)');
        % title(sprintf('%s Surrogate Model, Method: %s', stageName, surrogateType));
        % colorbar; 
        % grid on; 
        % view(45,25);
        % 
        % xlim([0 200]);
        % ylim([0 60])
        % zlim([0 1]);
        % 
        % boxColor = "r";
        % 
        % drawRefinementBox(refMin, refMax, boxColor);

        % Actual sample points in white
        nexttile;
        scatter3(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3), 60, 'w', 'filled', 'MarkerEdgeColor','k');
        hold on;
        
        xlabel('x1 (Advert)'); 
        ylabel('x2 (Target)'); 
        zlabel('x3 (Skew)');
        title(sprintf('%s Sample Space', stageName));
        grid on; 
        view(45,25);

        xlim([0 200]);
        ylim([0 60])
        zlim([0 1]);

        if (strcmpi(surrogateType, 'poly'))
            nexttile
        end

end