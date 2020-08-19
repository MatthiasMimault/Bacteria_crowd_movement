%% Plot script
% Run figure plot for the current loaded data with graphic parameters
% Last update 28/04/2020
addpath('..\Include')
addpath('..\Source')
addpath('..\Outputs')

%% Graphics settings
type = 'png'; % png, fig
scale = [10 200]; %min max bacteria density


%% Run
plotDensityDistribution(dataRoot, caseName, fileName ,nx, type, scale)