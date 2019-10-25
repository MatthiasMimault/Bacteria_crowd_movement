%% Executable script for Bacteria Movement model
% V3.0 - New structure
% Current: 
% Pending: obstacle generation
close all
clear

%% 0. Settings
% Names
caseName = 'DEV-BdyFlow';
caseType = 'A';
caseDate = '1025';
runnb = '2x';
suffix = 'HomogObstaDiff-T1000B01C001';
global R Dx Dy typeAt typeInit typeSrc typeObs typeVel A C epsilon BactValue

% Physical parameters
% T is time in s, A is bacteria average velocity, C is diffusion coeff.
% R is /// mm, epsilon is interaction strength
T = 1000; A = 0.03; C = 0.01; 
R = 2; epsilon = 0.5;
BactValue = 0.1;
typeAt = 'root'; %up, root
typeSrc = 'bottom';
typeObs = 'particles'; % none
typeVel = 'att-src'; % static, att-src
typeInit = 'square'; % none square

% Run parameters
debug = 'false'; % true false
Nfiles = 10; 

% Numerical parameters
nx = 400;
CFL = 0.9;
ny = nx;

% Dimensions
Domain = [-25, 25, -25, 25];
Space = [-20, 20, -20, 20];
Attractant = [-1 1 10 25]; %[-20 20 20 25] [-1 1 10 25]
InitBact = [-20, 20, -20, 20]; %[-20, 10, -10, 10]
Source = [-20 20 -22 -20]; 

% Folders
addpath('..\Data')
addpath('..\Include')
addpath('..\Figures')
addpath('..\Source')
if strcmp(debug, 'true') 
    addpath('..\Debug')
end

cd ..

if strcmp(debug, 'true') 
%     dataFolder = ['..\Debug\' caseDate '-' caseType runnb '-' num2str(nx)];
    dataFolder = ['\Debug\' caseName '\' caseDate '-' caseType runnb ...
        '-' num2str(nx) '-' suffix];
else
    dataFolder = ['\Data\' caseDate '-' caseType runnb '-' num2str(nx)];
end
figureFolder = ['\Figures\F' caseDate '-' caseType runnb '-' num2str(nx)];
pngFolder =    ['\Figures\P' caseDate '-' caseType runnb '-' num2str(nx)];

mkdir(dataFolder)
mkdir(figureFolder)
mkdir(pngFolder)

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
b = fDensityGeneration(X,Y,BactValue,InitBact);

% Generation convolution
[E,Ex,Ey] = fConvKernelGeneration();

% Generation velocity field
[Vxo, Vyo] = fVelocityGeneration(X, Y, domAt, domBd, domSrc);
Vx = Vxo;
Vy = Vyo;

% Graphics
Axis = Space;

% Save initial parameters
save([dataFolder '\' caseType runnb '-init'],'domDef','domBd',...
    'domAt', 'domSrc',...
    'Nt','T','Nfiles', 'Axis', 'caseName','dataFolder','nx')
save([dataFolder '\' caseType runnb '-000'],...
     'X','Y','b','Dx','Dy','Vx','Vy','tsave')

% Estimation time - start timer
tic
date = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'));
fprintf(['Start : ' date '\n']);

%% 2. Loop
for n = 1:Nt    
    %% Update
    % Velocity update
    [Vx,Vy] = fVelocityUpdate(b,Vxo,Vyo,Ex,Ey,domBd);
    
    % Density update
    b_final = fDensityUpdate(Dt,b,Vx,Vy,domBd,domSrc);
    
    % Total mass update
    B(n) = sum(sum(Dx*Dy*b));
    
    %% Save
    if n*Dt>=TT(itt+1)
        Dt_save = TT(itt+1)-Dt*(n-1);
        tsave =  Dt*(n-1)+Dt_save;
        
        b_save = fDensityUpdate(Dt_save,b,Vx,Vy,domBd,domSrc);
        plotSurf(X,Y,b_save,Vx,Vy,Axis)
        
        s = sprintf('%03s',num2str(itt,'%d'));
        save([dataFolder '\' caseType runnb '-' s],...
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
save([dataFolder '\' caseType runnb '-init'],...
    'B','Tcomp','tt','TT','-append')
