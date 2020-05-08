%% Plot function
% If released, remove pD and pD2
function plotDensityDistribution(dataRoot,caseName,fileName,nx,type,scale)
warning('off','MATLAB:contour:ConstantData')
close all
% Initialisation
% nameFolder = [nameFolder '\'];
load([dataRoot '\Data-' caseName '\' fileName '-init'],...
    'Axis', 'Nfiles','domBd','domSrc');

if strcmp(type,'png')
    folderPngs = [dataRoot '\Png-' caseName '-' fileName '-' num2str(nx)];
    if exist(folderPngs, 'dir') == 0
        mkdir(folderPngs)
    end
    
else
    folderFigs = [dataRoot '\Fig-' caseName '-' fileName '-' num2str(nx)];
    
    if exist(folderFigs, 'dir') == 0
        mkdir(folderFigs)
    end
end


% Time estimation
tic
date = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'));
fprintf(['Start : ' date '\n']);

for n = 0:Nfiles-1 
    % Load
    s = sprintf('%03s',num2str(n,'%d'));
    load([dataRoot '\Data-' caseName '\' fileName '-' s],...
        'X','Y','b','Vx','Vy')
    % Plot  
    plotSurf(X,Y,b,Vx,Vy,Axis, scale)     
    
    % Plot boundary
    hold on
    plotBd(X,Y,domBd+domSrc,Axis)
    hold off
    
    % Save
    if strcmp(type,'png')
        print([folderPngs '\' fileName '-' num2str(nx) '-' s],'-dpng')
    else
        savefig([folderFigs '\' fileName '-' num2str(nx) '-' s '.fig'])
    end
%     
    % print([nameFig 'Png_' s],'-dpng')
    
    % Estimation time
    plotTime(toc, n, Nfiles)
    pause(0.01)
end
plotFinalTime();