% sam test space
clear all, clc

filename = "ACTUALTEST.csv";

minAdvert = 0;
maxAdvert = 200;
minTarget = 0;
maxTarget = 60;
minSkew = 0;
maxSkew = 1;
maxVariables = [maxAdvert, maxTarget, maxSkew];
minVariables = [minAdvert, minTarget, minSkew];

colourRed = [1 0.1 0.1];
colourOrange = [1.0 0.5 0.1];
colourGreen = [0.1 1 0.1];

optsGA = optimoptions('ga', ...  
    'MaxGenerations', 1000, ... 
    'PopulationSize', 10000);  

dimensions = 3;
dstarIterations = 10000;
nodes = 30;
saveDStarDistribution = false;

% blackBoxFunc = @sphere3;
% blackBoxFunc = @griewank3;

% plotTrueExampleFunction(blackBoxFunc);


samples = [];
mcsi_values = [];

%figure initialize
% figAll = figure('Name','All Surrogates Comparison','NumberTitle','off');
% tiledlayout(figAll, 2, 2); %one Master popup for the plotting
% 
% % 1st Iteration.
% [pointsMatrix, bestDStar] = bestDStarPointDistribution(dimensions, 30, dstarIterations, saveDStarDistribution);
% tempSamplePoints = denormaliseMatrix(pointsMatrix, minVariables, maxVariables);
% % samples = [samples; tempSamplePoints];
% 
% % temp_mcsi_values = blackBoxFunc(samples); % or take interactive samples
% [samples, mcsi_values] = manualInputBlackbox(filename, tempSamplePoints);
% % %================================================================================
% % mcsi_values = [mcsi_values; temp_mcsi_values];
% 
% % 1st model
% zoomFactor = 1.25;
% [refMin, refMax] = fitSurrogateAndZoomArea(samples, mcsi_values, zoomFactor, true, "30 node Exploration", colourRed);


figAll2 = figure('Name','All Surrogates Comparison','NumberTitle','off');
tiledlayout(figAll2, 2, 2); %one Master popup for the plotting
% 2nd iteration

% normSamples = normaliseMatrix(samples, minVariables, maxVariables);
% normRefMin = normaliseMatrix(refMin, minVariables, maxVariables);
% normRefMax = normaliseMatrix(refMax, minVariables, maxVariables);
% [pointsMatrix, bestDStar, numExistingSamplesInRefArea] = refinementAreaPointDistribution(samples, refMin, refMax, minVariables, maxVariables, dimensions, 5, dstarIterations);
% % [pointsMatrix, bestDStar] = bestDStarPointDistribution(dimensions, 5, dstarIterations, saveDStarDistribution);
% tempSamplePoints = denormaliseMatrix(pointsMatrix, minVariables, maxVariables);
% % samples = [samples; tempSamplePoints];
% 
% % [samples, mcsi_values] = manualInputBlackbox(filename, tempSamplePoints);
% % temp_mcsi_values = blackBoxFunc(samples); % or take interactive samples
% [samples, mcsi_values] = manualInputBlackbox(filename, tempSamplePoints);
% mcsi_values = [mcsi_values; temp_mcsi_values];

% [samples, mcsi_values] = loadCSVSamples(filename);
% 
% [samples, mcsi_values] = loadCSVSamples(filename);
% zoomFactor = 0.9;
% % % [refMin, refMax] = fitSurrogateAndZoomArea(samples, mcsi_values, zoomFactor, true, "testStageName2");
% [refMin, refMax] = fitSurrogateAndZoomArea(samples, mcsi_values, zoomFactor, true, "Additional 5 Nodes in Refined Zone", colourRed)

% refMin = [ 50.4272   24.2835         0];
% refMax = [200.0000   60.0000    0.6545];
% [samples, mcsi_values] = loadCSVSamples(filename);
% 
% figAll3 = figure('Name','All Surrogates Comparison','NumberTitle','off');
% tiledlayout(figAll3, 2, 2); %one Master popup for the plotting
% % 3rd iteration
% % normSamples = normaliseMatrix(samples, minVariables, maxVariables);
% % normRefMin = normaliseMatrix(refMin, minVariables, maxVariables);
% % normRefMax = normaliseMatrix(refMax, minVariables, maxVariables);
% % [samples, mcsi_values] = loadCSVSamples(filename);
% [pointsMatrix, bestDStar, numExistingSamplesInRefArea] = refinementAreaPointDistribution(samples, refMin, refMax, minVariables, maxVariables, dimensions, 5, dstarIterations);
% % [pointsMatrix, bestDStar] = bestDStarPointDistribution(dimensions, 5, dstarIterations, saveDStarDistribution);
% tempSamplePoints = denormaliseMatrix(pointsMatrix, minVariables, maxVariables);
% % samples = [samples; tempSamplePoints];
% 
% % [samples, mcsi_values] = manualInputBlackbox(filename, tempSamplePoints);
% % temp_mcsi_values = blackBoxFunc(samples); % or take interactive samples
% [samples, mcsi_values] = manualInputBlackbox(filename, tempSamplePoints);
% % mcsi_values = [mcsi_values; temp_mcsi_values];
% 
% zoomFactor = 0.6;
% % [refMin, refMax] = fitSurrogateAndZoomArea(samples, mcsi_values, zoomFactor, true, "testStageName3");
% [refMin, refMax] = fitSurrogateAndZoomArea(samples, mcsi_values, zoomFactor, true, "Additional 5 Nodes in Further Refined Zone", colourOrange);


figAll4 = figure('Name','All Surrogates Comparison','NumberTitle','off');
tiledlayout(figAll4, 2, 2); %one Master popup for the plotting
% 3rd iteration
% normSamples = normaliseMatrix(samples, minVariables, maxVariables);
% normRefMin = normaliseMatrix(refMin, minVariables, maxVariables);
% normRefMax = normaliseMatrix(refMax, minVariables, maxVariables);
% [pointsMatrix, bestDStar, numExistingSamplesInRefArea] = refinementAreaPointDistribution(samples, refMin, refMax, minVariables, maxVariables, dimensions, 5, dstarIterations);
% [pointsMatrix, bestDStar] = bestDStarPointDistribution(dimensions, 5, dstarIterations, saveDStarDistribution);
% tempSamplePoints = denormaliseMatrix(pointsMatrix, minVariables, maxVariables);
% samples = [samples; tempSamplePoints];

% [samples, mcsi_values] = manualInputBlackbox(filename, tempSamplePoints);
% temp_mcsi_values = blackBoxFunc(samples); % or take interactive samples
% [samples, mcsi_values] = manualInputBlackbox(filename, tempSamplePoints);
% mcsi_values = [mcsi_values; temp_mcsi_values];


[samples, mcsi_values] = loadCSVSamples(filename);
zoomFactor = 0.4;
% [refMin, refMax] = fitSurrogateAndZoomArea(samples, mcsi_values, zoomFactor, true, "testStageName3");
[refMin, refMax] = fitSurrogateAndZoomArea(samples, mcsi_values, zoomFactor, true, "Additional 5 Nodes in Further Refined Zone",colourOrange);


figAll5 = figure('Name','All Surrogates Comparison','NumberTitle','off');
tiledlayout(figAll5, 2, 2); %one Master popup for the plotting
% 3rd iteration
% normSamples = normaliseMatrix(samples, minVariables, maxVariables);
% normRefMin = normaliseMatrix(refMin, minVariables, maxVariables);
% normRefMax = normaliseMatrix(refMax, minVariables, maxVariables);
[pointsMatrix, bestDStar, numExistingSamplesInRefArea] = refinementAreaPointDistribution(samples, refMin, refMax, minVariables, maxVariables, dimensions, 5, dstarIterations);
% [pointsMatrix, bestDStar] = bestDStarPointDistribution(dimensions, 5, dstarIterations, saveDStarDistribution);
tempSamplePoints = denormaliseMatrix(pointsMatrix, minVariables, maxVariables);
% samples = [samples; tempSamplePoints];

% [samples, mcsi_values] = manualInputBlackbox(filename, tempSamplePoints);
% temp_mcsi_values = blackBoxFunc(samples); % or take interactive samples
[samples, mcsi_values] = manualInputBlackbox(filename, tempSamplePoints);
% mcsi_values = [mcsi_values; temp_mcsi_values];

zoomFactor = 0.2;
% [refMin, refMax] = fitSurrogateAndZoomArea(samples, mcsi_values, zoomFactor, true, "testStageName3");
[refMin, refMax] = fitSurrogateAndZoomArea(samples, mcsi_values, zoomFactor, true, "Additional 5 Nodes in Further Refined Zone", colourGreen);