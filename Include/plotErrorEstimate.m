function [e, Dx2] = plotErrorEstimate(nameFolder,nameFile,nx1,nx2)
close all
%% Initialisation
% fast swap
if nx1>nx2
    nxt = nx1;
    nx1 = nx2;
    nx2 = nxt;
end

% Assume that the coarser discretisation is nx1, and identical final T and
% file number
% folderFigs1 = [nameFolder '\Figs-' nameFile '-' num2str(nx1,'%d')];

folderData1 = [nameFolder '\Data-' nameFile '-' num2str(nx1,'%d')];
load([folderData1 '\' nameFile '-' num2str(nx1,'%d') '-init'], 'Nt',...
    'dN', 'Axis');
% Nt1 = Nt;
% dN1 = dN;

% folderFigs2 = [nameFolder '\Figs-' nameFile '-' num2str(nx2,'%d')];
folderData2 = [nameFolder '\Data-' nameFile '-' num2str(nx2,'%d')];
load([folderData2 '\' nameFile '-' num2str(nx2,'%d') '-init'], 'Nt',...
    'dN', 'Nfiles');
% Nt2 = Nt;
% dN2 = dN;

% File index
nF = Nfiles;

% Time estimation
tic
date = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'));
fprintf(['Start : ' date '\n']);

errorL1 = size(1:nF+1);

for n = 0:nF-1  
    % Load b2
    s = sprintf('%03s',num2str(n,'%d'));
    load([folderData2 '\' nameFile '-' num2str(nx2,'%d') '-' s],...
        'X','Y','b','Dx','Dy');
    X2 = X; Y2 = Y; b2 = b; Dx2 = Dx; Dy2 = Dy;% TT2 = TT
%     figure(1)
%     contourf(X2,Y2,b2) 
%     colorbar
    
    % Load and upscale of b1
    s = sprintf('%03s',num2str(n,'%d'));
    load([folderData1 '\' nameFile '-' num2str(nx1,'%d') '-' s],...
        'X','Y','b');
    b1 = interp2(X,Y,b,X2,Y2,'nearest');% TT1 = TT
    b1(isnan(b1)) = 0; % NaN values on the boundary
%     figure(2)
%     contourf(X2,Y2,b1) 
%     colorbar
%     figure(4)
%     contourf(X,Y,b) 
%     colorbar
    
    % Error computation
    errorL1(n+1) = sum(sum(abs(Dx2*Dy2*b2-Dx2*Dy2*b1)));
    
    
%     % Plot  
% figure(3)
%     contourf(X2,Y2,b2-b1)    
%     colorbar
%     pause
%     % Save
%     savefig([folderFigs1 '\' nameFile '-' num2str(nx1,'%d') '-' s '.fig'])
    % print([nameFig 'Png_' s],'-dpng')
    
    % Estimation time
%     plotTime(toc, n, nF)
    pause(0.01)
end

% plot(T/nF:T/nF:T, errorL1)
% e = errorL1(end);
e = errorL1

plotFinalTime();