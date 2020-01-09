%% Executable script for Bacteria Movement model
% V3.0 - New structure, bacteria max density: 200 cells per mm^2
% V3.1 - New folder structure for figures: ready for sensitivity analysis
% Current: 
% Pending: Automated case generation, progress bar, plot to post
close all
clear

%% 0. Settings
% Names
caseTitle = 'SNS-woPtc';
caseDate = '0107';
caseType = 'B';
runnb = '6';
suffix = 'T500R05-e1C0001';
global T R Dx Dy nx CFL typeAt typeInit typeSrc typeObs typeVel A C epsilon...
    BactValue debug Ro ro

% Physical parameters
% T is time in s, A is bacteria average velocity, C is diffusion coeff.
% R is interaction radius in mm, epsilon is interaction strength
T = 500; A = 0.03; C = 0.001; 
R = 0.5; epsilon = 1;
BactValue = 20;
% particles - porosity 65%
Ro = 1.18; ro = 1;

typeAt = 'root'; %up, root
typeSrc = 'bottom';
typeObs = 'none'; % none, particles
typeVel = 'att-src'; % static, att-src, none, adv-src
typeInit = 'square'; % none square

% Run parameters
debug = 'false'; % true false
Nfiles = 10; 
typeRepulsion = 0;

% Numerical parameters
nx = 800;
CFL = 0.9;
ny = nx;

% Dimensions
Domain = [-4.25, 4.5, -2.2, 6];
Space = [-3.25, 3.5, -1.2, 5];
Attractant = [-0.2 0.2 4 6]; %[-20 20 20 25] [-1 1 10 25]
InitBact = [-3.25, 3.5, -1.2, 5]; %[-20, 10, -10, 10]
Source = [-3.25, 3.5 -2.2 -1.2]; 

% Directories
addpath('..\Include')
addpath('..\Source')
addpath('..\Cases')
addpath('..\Source')
% [dataRoot, dataFolder, figureFolder, pngFolder] = fFolderMaker( ...
%     caseName, caseDate, caseType, runnb, nx, suffix);
[dataRoot, caseName, fileName] = fFolderMaker2( ...
    caseTitle, caseDate, caseType, runnb, nx, suffix);

%% 1. Initialisation grid and dependant parameters
[X,Y,Dx,Dy] = fGridGeneration(nx,ny,Domain);

[domAt, domBd, domDef, domSrc] ...
    = fRegionGeneration(X, Y, Space, Attractant, Source);

Dt = min(min(Dx,Dy)/A,Dx*Dx/2/C)*CFL;
disp(['Time step is ' num2str(Dt)]);
tt = 0:Dt:T;
TT = 0:T/Nfiles:T;
tsave = 0;
Nt = length(tt)-1;
B = zeros(1,length(1:Nt));
itt = 1;

% Generation density
% >>> Define a check that InitBact is strictly included in Space
b = fDensityGeneration(X,Y,BactValue,InitBact).*domDef;

% Generation convolution
[E,Ex,Ey] = fConvKernelGeneration();

% Generation velocity field
[Vxo, Vyo] = fVelocityGeneration(X, Y, domAt, domBd, domSrc);
Vx = Vxo;
Vy = Vyo;

% Graphics
Axis = Space;

% Save initial parameters
save([dataRoot '\Data-' caseName '\' fileName '-init'],'domDef','domBd',...
    'domAt', 'domSrc',...
    'Nt','T','Nfiles', 'Axis','dataRoot', 'fileName', 'caseName','nx')
save([dataRoot '\Data-' caseName '\' fileName '-000'],...
     'X','Y','b','Dx','Dy','Vx','Vy','tsave')

% Estimation time - start timer
tic
date = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'));
fprintf(['Start : ' date '\n']);

%% 2. Loop
for n = 1:Nt    
    %% Update
    % Velocity update
%     if n*Dt>500
%         typeRepulsion = 1;
%     end
    [Vx,Vy] = fRootFieldUpdate(Vxo,Vyo,typeRepulsion);
    [Vx,Vy] = fVelocityUpdate(b,Vx,Vy,Ex,Ey,domBd);
%     [Vx,Vy] = fVelocityUpdate(b,Vxo,Vyo,Ex,Ey,domBd);
    
    % Density update
    b_final = fDensityUpdate(Dt,b,Vx,Vy,domBd,domSrc);
    
    % Total mass update
    B(n) = sum(sum(Dx*Dy*b));
    
    %% Save
    if n*Dt>=TT(itt+1)
        Dt_save = TT(itt+1)-Dt*(n-1);
        tsave =  Dt*(n-1)+Dt_save;
        
        b_save = fDensityUpdate(Dt_save,b,Vx,Vy,domBd,domSrc);
%         plotSurf(X,Y,b_save,Vx,Vy,Axis)
        
        s = sprintf('%03s',num2str(itt,'%d'));
        save([dataRoot '\Data-' caseName '\' fileName '-' s],...
            'X','Y','b','Dx','Dy','Vx','Vy','tsave')

        % Clean-up
        plotTime(toc, n, Nt)
        pause(0.01)
        itt = itt + 1;
    end
    
    %% Final
    b = b_final;

end

%% 3. Post-traitment
plotFinalTime();
Tcomp = toc;
tt = Dt*(1:Nt);
save([dataRoot '\Data-' caseName '\' fileName '-init'],...
    'B','Tcomp','tt','TT','-append')
postLogGeneration([dataRoot '\Data-' caseName '\'], fileName,...
    Domain, Space, Attractant, Source, InitBact)

