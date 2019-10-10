%% Executable script for Bacteria Movement model
% V3.0 - New structure
close all
clear

%% 0. Settings
% Names
caseName = 'DEV-BdyFlow';
caseType = 'A';
caseDate = '1010';
runnb = '1';
suffix = 'SimpleAdvection';

% Physical parameters

% Run parameters
debug = 'true';

% Numerical parameters
nx = 100;

% Dimensions

% Folders
addpath('..\Data')
addpath('..\Include')
addpath('..\Figures')
addpath('..\Source')
if strcmp(debug, 'true') 
    addpath('..\Debug')
end

if strcmp(debug, 'true') 
    dataFolder = ['..\Debug\' caseDate '-' caseType runnb '-' num2str(nx)];
else
    dataFolder = ['..\Data\' caseDate '-' caseType runnb '-' num2str(nx)];
end
figureFolder = ['..\Figures\F' caseDate '-' caseType runnb '-' num2str(nx)];
pngFolder =    ['..\Figures\P' caseDate '-' caseType runnb '-' num2str(nx)];

mkdir(dataFolder)
mkdir(figureFolder)
mkdir(pngFolder)
