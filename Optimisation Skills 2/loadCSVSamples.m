% our format of storing data
function [denormalisedPointsMatrix, MCSI_values] = loadCSVSamples(path)
    csvFile = readtable(path);
    
    x1 = csvFile.x1;
    x2 = csvFile.x2;
    x3 = csvFile.x3;
    MCSI_values = csvFile.MCSI;
    denormalisedPointsMatrix = [x1, x2, x3];
end