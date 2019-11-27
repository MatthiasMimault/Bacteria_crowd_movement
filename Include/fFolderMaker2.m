function [dataRoot, caseName, fileName] = fFolderMaker2(...
    caseTitle, caseDate, caseType, runnb, nx, suffix)
global debug

mkdir('..\Cases')
if strcmp(debug, 'true') 
    mkdir('..\Debug')
end

cd ..

if strcmp(debug, 'true') 
%     dataFolder = ['..\Debug\' caseDate '-' caseType runnb '-' num2str(nx)];
    dataRoot = ['Debug\' caseTitle];
else
    dataRoot = ['Cases\' caseTitle];
end
mkdir(dataRoot)

caseName = [caseDate '-' caseType runnb '-' num2str(nx) '-' suffix];
fileName = [caseType runnb];
dataFolder = [dataRoot '\Data-' caseName];
mkdir(dataFolder)