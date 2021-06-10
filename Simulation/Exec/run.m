%% Executable script for Bacteria Crowd Movement model
%     Copyright (C) 2021  Matthias Mimault
%     Contact matthias.mimault -at- hutton.ac.uk
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
close all
clear

%% 0. Settings
% Names
caseTitle = 'Obstacles';
caseDate = '0906'; %MD
suffix = 'simple_1';
global T Dx Dy nx CFL A D Ro R BactValue domAt domBd domSrc domDef

% Load default parameters
input_parameters

% Custom parameters
Nfiles = 5; 
D = 0.05;
R = 1;
nx = 200;

typeObs = 'particles'; % none, particles
Ro = 4;
Domain = [-10, 10, -14.5, 6];
Space = [-8, 8, -13.2, 5];
InitBact = Space;

%% Nothing to be modified below
typeRepulsion = 0;

% Directories
addpath('..\Include')
addpath('..\Source')
addpath('..\Outputs')
[dataRoot, caseName, fileName] = fFolderMaker35( ...
    caseTitle, caseDate, nx, suffix);

% Licence disclaimer
disp('Bacteria crowd mvt. model  Copyright (C) 2021  Matthias Mimault')
disp('This program comes with ABSOLUTELY NO WARRANTY. This is free') 
disp('software and you are welcome to redistribute it under certain')
disp('conditions.')

%% 1. Initialisation grid and dependant parameters
[X,Y,Dx,Dy] = fGridGeneration(nx,Domain);

[domAt, domBd, domDef, domSrc] ...
    = fRegionGeneration(X, Y, Space, Attractant, Source);

% Dt = min(min(Dx,Dy)/A/(1+epsilon),Dx*Dx/2/C)*CFL;
Dt = min(Dx/A,Dx*Dx/2/D)*CFL;
disp(['Simulation name is ' caseTitle '/' caseName])
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
postLogGeneration([dataRoot '\Data-' caseName '\'], caseName,...
    Domain, Space, Attractant, Source, InitBact)

plot_gener

