% sam test space
clear all, clc

filename = "sam.csv";

minAdvert = 0;
maxAdvert = 200;
minTarget = 0;
maxTarget = 60;
minSkew = 0;
maxSkew = 1;
maxVariables = [maxAdvert, maxTarget, maxSkew];
minVariables = [minAdvert, minTarget, minSkew];

optsGA = optimoptions('ga', ...  
    'MaxGenerations', 100, ... 
    'PopulationSize', 1000);  

dimensions = 3;
dstarIterations = 5000;
nodes = 30;
saveDStarDistribution = false;

% blackBoxFunc = @sphere3;
blackBoxFunc = @griewank3;

plotTrueExampleFunction(blackBoxFunc);

samples = [];
mcsi_values = [];


% 1st Iteration
[pointsMatrix, bestDStar] = bestDStarPointDistribution(dimensions, 50, dstarIterations, saveDStarDistribution);
tempSamplePoints = denormaliseMatrix(pointsMatrix, minVariables, maxVariables);
samples = [samples; tempSamplePoints];

temp_mcsi_values = blackBoxFunc(samples); % or take interactive samples
% [samples, mcsi_values] = manualInputBlackbox(filename, tempSamplePoints);
% %================================================================================
mcsi_values = [mcsi_values; temp_mcsi_values];

% 1st model
zoomFactor = 0.9;
[refMin, refMax] = fitSurrogateAndZoomArea(samples, mcsi_values, zoomFactor, true, "testStageName");



% 2nd iteration
[pointsMatrix, bestDStar] = bestDStarPointDistribution(dimensions, 5, dstarIterations, saveDStarDistribution);
tempSamplePoints = denormaliseMatrix(pointsMatrix, refMin, refMax);
samples = [samples; tempSamplePoints];

% [samples, mcsi_values] = manualInputBlackbox(filename, tempSamplePoints);
temp_mcsi_values = blackBoxFunc(tempSamplePoints); % or take interactive samples
mcsi_values = [mcsi_values; temp_mcsi_values];

zoomFactor = 0.6;
% [refMin, refMax] = fitSurrogateAndZoomArea(samples, mcsi_values, zoomFactor, true, "testStageName2");
[refMin, refMax] = fitSurrogateAndZoomArea(tempSamplePoints, temp_mcsi_values, zoomFactor, true, "testStageName2");



% 3rd iteration
[pointsMatrix, bestDStar] = bestDStarPointDistribution(dimensions, 5, dstarIterations, saveDStarDistribution);
tempSamplePoints = denormaliseMatrix(pointsMatrix, refMin, refMax);
samples = [samples; tempSamplePoints];

% [samples, mcsi_values] = manualInputBlackbox(filename, tempSamplePoints);
temp_mcsi_values = blackBoxFunc(tempSamplePoints); % or take interactive samples
mcsi_values = [mcsi_values; temp_mcsi_values];

zoomFactor = 0.4;
% [refMin, refMax] = fitSurrogateAndZoomArea(samples, mcsi_values, zoomFactor, true, "testStageName3");
[refMin, refMax] = fitSurrogateAndZoomArea(tempSamplePoints, temp_mcsi_values, zoomFactor, true, "testStageName3");