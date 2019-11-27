function [dataRoot, dataFolder, figureFolder, pngFolder] = fFolderMaker(...
    caseName, caseDate, caseType, runnb, nx, suffix)
global debug
addpath('..\Data')
addpath('..\Figures')
if strcmp(debug, 'true') 
    addpath('..\Debug')
end

cd ..

if strcmp(debug, 'true') 
%     dataFolder = ['..\Debug\' caseDate '-' caseType runnb '-' num2str(nx)];
    dataRoot = ['Debug\' caseName];
else
    dataRoot = ['Data\' caseName];
end
mkdir(dataRoot)

if strcmp(debug, 'true') 
%     dataFolder = ['..\Debug\' caseDate '-' caseType runnb '-' num2str(nx)];
    dataFolder = ['Debug\' caseName '\' caseDate '-' caseType runnb ...
        '-' num2str(nx) '-' suffix];
else
    dataFolder = ['Data\' caseName '\' caseDate '-' caseType runnb ...
        '-' num2str(nx) '-' suffix];
end
figureFolder = ['Figures\F' caseDate '-' caseType runnb '-' num2str(nx)];
pngFolder =    ['Figures\P' caseDate '-' caseType runnb '-' num2str(nx)];

mkdir(dataFolder)
mkdir(figureFolder)
mkdir(pngFolder)