function [dataRoot, caseName, fileName] = fFolderMaker35(...
    caseTitle, caseDate, nx, suffix)
global dev

if exist('..\Outputs', 'dir')==0
    mkdir('..\Outputs')
end

switch dev
    case 'stale'
        dataRoot = ['..\Outputs\Stl-' caseTitle];
    case 'dev'
        dataRoot = ['..\Outputs\Dev-' caseTitle];
    otherwise
        dataRoot = ['..\Outputs\' caseTitle];
end

if exist(dataRoot, 'dir') == 0
    mkdir(dataRoot)
end


caseName = [caseDate '-' suffix '-' num2str(nx)];
fileName = suffix;
dataFolder = [dataRoot '\Data-' caseName];
if exist(dataFolder, 'dir') == 0
    mkdir(dataFolder)
end