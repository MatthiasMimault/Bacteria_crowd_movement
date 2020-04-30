function [dataRoot, caseName, fileName] = fFolderMaker(...
    caseTitle, caseDate, runnb, nx, suffix)
global dev

switch dev
    case {'stale', 'dev'}
        if exist('..\Dev', 'dir') == 0
            mkdir('..\Dev')
        end
    otherwise
        if exist('..\Cases', 'dir') == 0
            mkdir('..\Cases')
        end
end

switch dev
    case 'stale'
        dataRoot = ['..\Dev\Stl-' caseTitle];
    case 'dev'
        dataRoot = ['..\Dev\Dev-' caseTitle];
    otherwise
        dataRoot = ['..\Cases\' caseTitle];
end

if exist(dataRoot, 'dir') == 0
    mkdir(dataRoot)
end


caseName = [caseDate '-' runnb '-' num2str(nx) '-' suffix];
fileName = [runnb];
dataFolder = [dataRoot '\Data-' caseName];
if exist(dataFolder, 'dir') == 0
    mkdir(dataFolder)
end