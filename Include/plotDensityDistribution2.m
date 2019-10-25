%% Plot function
function plotDensityDistribution2(nameFolder,nameFile,nx,type)
close all
% Initialisation
% nameFolder = [nameFolder '\'];
load([nameFolder '\' nameFile '-init'],...
    'Axis', 'Nfiles','domBd');

if strcmp(type,'png')
    folderPngs = [nameFolder '\Pngs-' nameFile '-' num2str(nx)];
    mkdir(folderPngs)
else
    folderFigs = [nameFolder '\Figs-' nameFile '-' num2str(nx)];
    mkdir(folderFigs)
end


% Time estimation
tic
date = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'));
fprintf(['Start : ' date '\n']);

for n = 0:Nfiles-1 
    % Load
    s = sprintf('%03s',num2str(n,'%d'));
    load([nameFolder '\' nameFile '-' s],...
        'X','Y','b','Vx','Vy')
    % Plot  
    plotSurf(X,Y,b,Vx,Vy,Axis)     
    
    % Plot boundary
    hold on
    plotBd(X,Y,domBd,Axis)
    hold off
    
    % Save
    if strcmp(type,'png')
        print([folderPngs '\' nameFile '-' num2str(nx) '-' s],'-dpng')
    else
        savefig([folderFigs '\' nameFile '-' num2str(nx) '-' s '.fig'])
    end
%     
    % print([nameFig 'Png_' s],'-dpng')
    
    % Estimation time
    plotTime(toc, n, Nfiles)
    pause(0.01)
end
plotFinalTime();