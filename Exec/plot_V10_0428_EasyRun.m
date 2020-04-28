%% Plot script
% Run figure plot for the current loaded data with graphic parameters
% Last update 28/04/2020

%% Graphics settings
type = 'png'; % png, fig
scale = [15 200]; %min max bacteria density


%% Run
plotDensityDistribution2(dataRoot, caseName, fileName ,nx, type, scale)