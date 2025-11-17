function output = manualBlackbox(samples,stage)
    
    [denormalisedPointsMatrix, MCSI_values] = loadCSVSamples("manualEntry.csv");

    switch stage

    case 'initialize'   % Haven't made any points yet, 30 samples is incoming
        MCSI_values = nan(size(samples,1),1); %Populates the field with NaN so its writable
        dataMatrix = [samples, MCSI_values];   % N×4
        T = array2table(dataMatrix,'VariableNames', {'x1','x2','x3','MCSI'});
        writetable(T, "manualEntry.csv"); %Write data to .csv
        
        while true %wait for the user to input MCSI values
            userInput = input('Enter MCSI values or type "done": ', 's');
            
            if strcmpi(userInput, 'done')
                break;   % Exit loop only when user is done inputting MCSI values
            else
                pause(1);    % Sleep for 1 second
            end
        end

        [~, MCSI_values] = loadCSVSamples("manualEntry.csv"); %Loads the users freshly inputted values
        output = MCSI_values; %Returns the users inputted MCSI
        return;

    case 'poly'   % 5 samples for poly refinement incoming
        allPoints   = [denormalisedPointsMatrix; samples];   %Appends new data with existing data   
        newMCSI     = [MCSI_values; nan(size(samples,1),1)];    %Appends new data with existing data

        dataMatrix = [allPoints, newMCSI];   
        T = array2table(dataMatrix,'VariableNames', {'x1','x2','x3','MCSI'});
        writetable(T, "manualEntry.csv");
        while true %wait for the user to input MCSI values
                userInput = input('Enter MCSI values then type "done": ', 's');
                
                if strcmpi(userInput, 'done')
                    break;   % Exit loop
                else
                    pause(1);    
                end
            end
    
            [~, MCSI_values] = loadCSVSamples("manualEntry.csv");
            output = MCSI_values;
            return;

    case 'poly2'   % 5 samples for poly refinement 2 incoming
        allPoints   = [denormalisedPointsMatrix; samples];      
        newMCSI     = [MCSI_values; nan(size(samples,1),1)];    

        dataMatrix = [allPoints, newMCSI];   
        T = array2table(dataMatrix,'VariableNames', {'x1','x2','x3','MCSI'});
        writetable(T, "manualEntry.csv");
        while true %wait for the user to input MCSI values
                userInput = input('Enter MCSI values or type "done": ', 's');
                
                if strcmpi(userInput, 'done')
                    break;   % Exit loop
                else
                    pause(1);    
                end
            end
    
            [~, MCSI_values] = loadCSVSamples("manualEntry.csv");
            output = MCSI_values;
            return;

    case 'rbf'  % 5 samples for RBF refinement incoming
        allPoints   = [denormalisedPointsMatrix; samples];     
        newMCSI     = [MCSI_values; nan(size(samples,1),1)];   

        dataMatrix = [allPoints, newMCSI];   
        T = array2table(dataMatrix,'VariableNames', {'x1','x2','x3','MCSI'});
        writetable(T, "manualEntry.csv");
        while true %wait for the user to input MCSI values
                userInput = input('Enter MCSI values or type "done": ', 's');
                
                if strcmpi(userInput, 'done')
                    break;   % Exit loop
                else
                    pause(1);   
                end
            end
    
            [~, MCSI_values] = loadCSVSamples("manualEntry.csv");
            output = MCSI_values;
            return;

    case 'rbf2'  % 5 samples for rbf refinement 2 incoming
        allPoints   = [denormalisedPointsMatrix; samples];      
        newMCSI     = [MCSI_values; nan(size(samples,1),1)];    

        dataMatrix = [allPoints, newMCSI];        
        T = array2table(dataMatrix,'VariableNames', {'x1','x2','x3','MCSI'});
        writetable(T, "manualEntry.csv");
        while true %wait for the user to input MCSI values
                userInput = input('Enter MCSI values or type "done": ', 's');
                
                if strcmpi(userInput, 'done')
                    break;   % Exit loop
                else
                    pause(1);    
                end
            end
    
            [~, MCSI_values] = loadCSVSamples("manualEntry.csv");
            output = MCSI_values;
            return;

        otherwise
        error('Unexpected number of MCSI sample values.')
end


end