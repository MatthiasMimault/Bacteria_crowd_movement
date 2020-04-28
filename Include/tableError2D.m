%% Plot function
% If released, remove pD and pD2
function tableError2D(dataRoot,caseName,fileName,nx,type)
warning('off','MATLAB:contour:ConstantData')
close all
% Initialisation
% nameFolder = [nameFolder '\'];
load([dataRoot '\Data-' caseName '\' fileName '-init'],...
    'Nfiles','A','C','Dx','Dt','TT');

%% Exact solution
be = @(x,y) exp(-A/C*sqrt(x.^2+y.^2));


% Time estimation
tic
date = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'));
fprintf(['Start : ' date '\n']);
err = zeros(1,Nfiles);

for n = 0:Nfiles-1
% for n = 0:5
    % Load
    s = sprintf('%03s',num2str(n,'%d'));
    load([dataRoot '\Data-' caseName '\' fileName '-' s],...
        'X','Y','b')
    
    err(n+1) = Dx*Dx*sum(sum(abs(be(X,Y)-b)));
    
    % Estimation time
    plotTime(toc, n, Nfiles)
    pause(0.01)
end
% plot(X(end/2,:),be(X(end/2,:)),'r--')
% figure
err(end)
plot(TT(1:end-1),err)
plotFinalTime();