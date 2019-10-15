%% Executable script for Bacteria Movement model
% V3.0 - New structure
% Pending: convolution grid generation
close all
clear

%% 0. Settings
% Names
caseName = 'DEV-BdyFlow';
caseType = 'A';
caseDate = '1010';
runnb = '1';
suffix = 'SimpleAdvection';
global Dx Dy typeAt typeSrc typeObs typeVel

% Physical parameters
% T is time in s, A is bacteria average velocity, C is diffusion coeff.
% R is /// mm, epsilon is interaction strength
T = 200; A = 0.03; C = 0; 
R = 2; epsilon = 0.5;
typeAt = 'up';
typeSrc = 'bottom';
typeObs = 'none';
typeVel = 'att-adv';

% Run parameters
debug = 'true';
Nfiles = 10; 

% Numerical parameters
nx = 400;
CFL = 0.9;
ny = nx;

% Dimensions
Domain = [-25, 25, -25, 25];
Space = [-20, 20, -20, 20];
Attractant = [-20 20 20 25];
InitBact = [-19.25, 19.25, -19.25, -5.25]; 
Source = [-20 20 -25 -20]; 

% Folders
addpath('..\Data')
addpath('..\Include')
addpath('..\Figures')
addpath('..\Source')
if strcmp(debug, 'true') 
    addpath('..\Debug')
end

if strcmp(debug, 'true') 
%     dataFolder = ['..\Debug\' caseDate '-' caseType runnb '-' num2str(nx)];
    dataFolder = ['..\Debug\' caseName '\' caseDate '-' caseType runnb ...
        '-' num2str(nx) '-' suffix];
else
    dataFolder = ['..\Data\' caseDate '-' caseType runnb '-' num2str(nx)];
end
figureFolder = ['..\Figures\F' caseDate '-' caseType runnb '-' num2str(nx)];
pngFolder =    ['..\Figures\P' caseDate '-' caseType runnb '-' num2str(nx)];

mkdir(dataFolder)
mkdir(figureFolder)
mkdir(pngFolder)

%% 1. Initialisation grid and dependant parameters
[X,Y,Dx,Dy] = fGridGeneration(nx,ny,Domain);

[domAt, domBd, domDef, domSrc] ...
    = fRegionGeneration(X, Y, Space, Attractant, Source);
% contourf(domAt+2*domBd+3*domSrc+7*domDef)
% colorbar

Dt = min(min(Dx,Dy)/A,Dx*Dx/2/C)*CFL
tt = 0:Dt:T;
TT = 0:T/Nfiles:T;
Nt = length(tt)-1;
B = zeros(1,length(1:Nt));
itt = 1;

% Convolution
% >>> Required generation convolution grid

% Velocity field generation
[Vxo, Vyo] = fVelocityGeneration(X, Y, domAt, domBd, domDef);
Vx = Vxo;
Vy = Vyo;

% Graphics
Axis = Space;

% Save initial parameters
save([dataFolder '\' caseType runnb '-init'],'domDef','domBd',...
    'domAt', 'domSrc',...
    'Nt','T','Nfiles', 'Axis', 'caseName','dataFolder','nx')
save([dataFolder '\' caseType runnb '-000'])

% Estimation time - start timer
tic
date = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'));
fprintf(['Start : ' date '\n']);

%% 2. Loop
% >>> Required loop

%% 3. Post-traitment
plotFinalTime();
Tcomp = toc;
tt = Dt*(1:Nt);
save([dataFolder '\' caseType runnb '-init'],...
    'B','Tcomp','tt','TT','-append')
