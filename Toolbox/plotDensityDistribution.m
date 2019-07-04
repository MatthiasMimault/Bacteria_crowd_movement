%% Plot function
function plotDensityDistribution(nameFolder,nameFile,nx,type)
close all
% Initialisation
folderData = [nameFolder '\Data-' nameFile];
load([folderData '\' nameFile '-init'],...
    'Axis', 'Nfiles');

if strcmp(type,'png')
    folderPngs = [nameFolder '\Pngs-' nameFile];
    mkdir(folderPngs)
else
    folderFigs = [nameFolder '\Figs-' nameFile];
    mkdir(folderFigs)
end


% Time estimation
tic
date = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'));
fprintf(['Start : ' date '\n']);

for n = 1:Nfiles 
    % Load
    s = sprintf('%03s',num2str(n,'%d'));
    load([folderData '\' nameFile '-' s],...
        'X','Y','b','Vx','Vy')
    % Plot  
    plotSurf(X,Y,b,Vx,Vy,Axis)      
    
    % Save
    if strcmp(type,'png')
        print([folderPngs '\' nameFile '-' s],'-dpng')
    else
        savefig([folderFigs '\' nameFile '-' s '.fig'])
    end
%     
    % print([nameFig 'Png_' s],'-dpng')
    
    % Estimation time
    plotTime(toc, n, Nfiles)
    pause(0.01)
end
plotFinalTime();