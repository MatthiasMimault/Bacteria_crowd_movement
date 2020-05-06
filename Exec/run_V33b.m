%% Executable script for Bacteria Movement model
% V3.0 - New structure, bacteria max density: 200 cells per mm^2
% V3.1 - New folder structure for figures: ready for sensitivity analysis
% V3.2 - New density correction: I original (non-normalised)
% V3.3 - New density units: 1e6 cells per mL, fixed numerical scheme, clean
% Current: 
% Pending: Automated case generation, progress bar, plot to post
close all
clear

%% 0. Settings
% Names
caseTitle = 'RUN-SimpleCase';
caseDate = '0506';
runnb = '6';
suffix = 'T5-k005D001';
global T Dx Dy nx CFL A D BactValue 

% input file
input_V33

%% Nothing to be modified below
typeRepulsion = 0;

% Directories
addpath('..\Include')
addpath('..\Source')
addpath('..\Cases')
[dataRoot, caseName, fileName] = fFolderMaker( ...
    caseTitle, caseDate, runnb, nx, suffix);

%% 1. Initialisation grid and dependant parameters
[X,Y,Dx,Dy] = fGridGeneration(nx,Domain);

[domAt, domBd, domDef, domSrc] ...
    = fRegionGeneration(X, Y, Space, Attractant, Source);

% Dt = min(min(Dx,Dy)/A/(1+epsilon),Dx*Dx/2/C)*CFL;
Dt = min(Dx/A,Dx*Dx/2/D)*CFL;
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
% contourf(X,Y,b); colorbar

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
    'domAt', 'domSrc','A','D','Dx','Dt',...
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
    [Vx,Vy] = fRootFieldUpdate(Vxo,Vyo,typeRepulsion);
    [Vx,Vy] = fVelocityUpdate(b,Vx,Vy,Ex,Ey,domBd);
    
    % Density update
    b_final = fDensityUpdate(Dt,b,Vx,Vy,domBd,domSrc);
    
    % Total mass update
    B(n) = sum(sum(Dx*Dy*b));
    
    %% Save
    if n*Dt>=TT(itt+1)
        Dt_save = TT(itt+1)-Dt*(n-1);
        tsave =  Dt*(n-1)+Dt_save;
        
        b_save = fDensityUpdate(Dt_save,b,Vx,Vy,domBd,domSrc);
        
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

