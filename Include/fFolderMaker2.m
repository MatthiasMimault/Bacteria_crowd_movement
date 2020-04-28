function [dataRoot, caseName, fileName] = fFolderMaker2(...
    caseTitle, caseDate, caseType, runnb, nx, suffix)
global debug

if exist('..\Cases', 'dir') == 0
    mkdir('..\Cases')
end
if strcmp(debug, 'true')     
    if exist('..\Debug', 'dir') == 0
        mkdir('..\Debug')
    end
end


if strcmp(debug, 'true') 
    dataRoot = ['..\Debug\' caseTitle];
else
    dataRoot = ['..\Cases\' caseTitle];
end
if exist(dataRoot, 'dir') == 0
    mkdir(dataRoot)
end


caseName = [caseDate '-' caseType runnb '-' num2str(nx) '-' suffix];
fileName = [caseType runnb];
dataFolder = [dataRoot '\Data-' caseName];
if exist(dataFolder, 'dir') == 0
    mkdir(dataFolder)
end