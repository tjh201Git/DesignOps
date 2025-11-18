% manual blockbox input
function [samplePoints, MCSI] = manualInputBlackbox(filename, samples)

    if isfile(filename)
        disp('File exists.');
    else
        disp('File does not exist, creating file.');
        fid = fopen(filename, 'w'); % Open the file for writing
        fprintf(fid, 'x1,x2,x3,MCSI');
        fclose(fid); % Close the file immediately after creating it
    end

    [existingSamples, existingMCSIs] = loadCSVSamples(filename);
    existingSamplesLen = length(existingSamples(:, 1));

    x1 = samples(:, 1);
    x2 = samples(:, 2);
    x3 = samples(:, 3);

    samplesLen = length(x1);

    MCSI_values = nan(samplesLen,1); %Populates the field with NaN so its writable

    fprintf("%d samples already exist. Asking for %d more samples\n", existingSamplesLen, samplesLen);

    for i = 1: samplesLen
        x1_ = x1(i);
        x2_ = x2(i);
        x3_ = x3(i);

        actualSampleNum = i + existingSamplesLen;
        fprintf("Sample %d (%.4f, %.4f, %.4f). ", actualSampleNum, x1_, x2_, x3_);

        msci = NaN;
        while isnan(msci)
            userInput = input('Enter MCSI value: ', 's');
            msci = str2double(userInput);
            if isnan(msci)
                disp('Invalid input. Please enter a numeric value.');
            end
        end

        % Invert it so that we get a good minimum
        msci = msci * -1;

        MCSI_values(i) = msci;

        row = [x1_, x2_, x3_, msci];
        
        % append the values
        writematrix(row, filename, 'WriteMode', 'append');
        fprintf("Appended %.4f %.4f %.4f %.4f to %s\n\n", row, filename);


    end
           
    
    samplePoints = [existingSamples; samples];
    MCSI   = [existingMCSIs; MCSI_values];
end