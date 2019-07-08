%% Main file for Bacterial movement model
% V1.70 - Attractive root, repulsive region, weno3
% V1.71 - Tentative automatic file naming aborted, correction of NsaveVect
% generation (validated mostly, no wrong case found for CFL 0.1-0.8)
% V1.72 - Correction of behaviour at boundary: Solved, but Velocity
% orientation troubling
close all
clear

%% Parameters
% Settings: Time T, Advection velocity A, Attractant trail B
% Obstacles: ro is the size (side) of a square obstacle and
% Ro the spacing between the center of two obstacles
% R is the radius of interaction for the perceived density
T = 100; A = 0.03; B = 1; 
R = 0.2;
Nfiles = 20; 
nameFolder = 'DBG-BoundaryStability-0407';
suffix = 'Advection';
runnb = '0';
% A standstill, B advection, C adv-obstacles, D attraction, E att-obstacles
% F Root alone
type = 'F'; 

% Numerical parameters
nx = 300;
CFL = 0.9;
ny = nx;

% nameFile = fCaseNaming(type,suffix,nx,nameFolder);
nameFile = [type runnb '-' suffix '-' num2str(nx)];

% Domain
Domain = [-25, 25, -25, 25];
Space = [-20, 20, -20, 20];

% Bacteria
InitSpace = [-15, 15, -15, 0];
InitValue = 0.5;

switch type
    case 'A'
        typeAt = 'none';
        typeObs = 'none';
        typeVel = 'stand';
    case 'B'
        typeAt = 'up';
        typeObs = 'none';
        typeVel = 'adv';
    case 'C'
        typeAt = 'up';
        typeObs = 'particles';
        typeVel = 'adv';
    case 'D'
        typeAt = 'up';
        typeObs = 'none';
        typeVel = 'att';
    case 'E'
        typeAt = 'up';
        typeObs = 'particles';
        typeVel = 'att';
    case 'F'
        typeAt = 'root';
        typeObs = 'none';
        typeVel = 'adv';
        
    otherwise
        disp('Unknown case. Abort')
        return
end


% Serial obstacles: Ro obstacle center spacing, ro obstacle size, re
% repelent region around obstacles
% V1.7 - Deprecated, moved to fBoundaryGeneration
% Ro = 0.5; ro = 0.03; re = 0.06;
% Obstacles = fObstacleGeneration(Space,Ro,ro);
% Obstacles = [];

% Graphics
Axis = Space;

%% Initialisation
% make directory
mkdir('Results')
nameFolder = ['Results\' nameFolder];
folderData = [nameFolder '\Data-' nameFile];
mkdir(nameFolder)
mkdir(folderData)
addpath('..\Toolbox')  

% Define regions
% Domain generation (2D) (Finite Volumes Cartesian grid)
% V2.0 Generation ok, to be tested. To move Dxy inside function
% V3.0 Breakdown of region creation with string options
% [X,Y,Dx,Dy,PhiBd,PhiAt,PhiDef] = ...
%     fGridGeneration(nx,ny, Domain, Space, Attractant, Obstacles);

[X,Y,Dx,Dy] = fGridGeneration(nx,ny,Domain);
PhiAt = fAttractantGeneration(X,Y,Attractant,typeAt);
PhiBd = fBoundaryGeneration(X,Y,PhiAt,Space,typeAt,typeObs);
PhiDef = ones(size(X))-PhiAt-PhiBd;

% Constants domain
Dt = min(Dx,Dy)*CFL;
tt = 0:Dt:T;
TT = 0:T/Nfiles:T;
Nt = length(tt);


%% Convolution
% bacteria interaction kernel
[E,Ex,Ey] = fConvolutionKernels(R,Dx,Dy);

%% Variables; Density b, Attractant p
b = fInitDatGeneration(InitSpace, InitValue, X, Y, PhiDef);
% b = 0.5*PhiDef; 
p = zeros(size(X));
Upx = zeros(size(X)); Umx = zeros(size(X));
Upy = zeros(size(X)); Umy = zeros(size(X));
    

%% Eikonal equation for P
% 10/06 - V1.22 Global equation, remove cropping
% 28/06 - V1.62 Zeros vector field for rest
% 03/07 - V1.72 Switch cases

switch typeVel
    case 'stand'
        Vxo = zeros(size(X));
        Vyo = zeros(size(X));
    case {'adv','att'}
        PhiC = fEikonalCost(X,Y,ones(size(X)),PhiBd,PhiAt);
        [Vxo,Vyo] = fDiffFlex(PhiC,PhiDef,Dx,Dy);
        Vn = sqrt(Vxo.^2+Vyo.^2+PhiBd+PhiAt);
        Vxo = -Vxo./Vn;
        Vxo(isnan(Vxo)) = 0;
        Vyo = -Vyo./Vn;
        Vyo(isnan(Vyo)) = 0;
        Vxoo = Vxo; Vyoo = Vyo;
    otherwise 
        Vxo = zeros(size(X));
        Vyo = zeros(size(X));
end


%% Estimation time
tic
date = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'));
fprintf(['Start : ' date '\n']);

%% Mass evolution
B = zeros(1,length(1:Nt));

%% Save initial parameters
nSaveVec = zeros(size(TT));
for i=2:length(TT)-1
   nFind = find(tt-TT(i)>=0 & abs(tt-TT(i))<Dt);
   nSaveVec(i) = nFind(1);
end
nSaveVec(end) = Nt;
NSave = 1;
dN = ceil(Nt/Nfiles);
save([folderData '\' nameFile '-init'],...
    'Nt','dN','T','Nfiles', 'Axis', 'nameFolder','nameFile','nx')
save([folderData '\' nameFile '-000'])

%% Loop
for n = 2:Nt
    % Half step x
    switch typeVel
        case 'att'
            Vx = fDirectionX(b,Vxo,Vyo,Ex,Ey,PhiDef); 
        otherwise
            Vx = Vxo;
    end
    btemp = fUpdateX(b,CFL,Vx, PhiDef, PhiAt, PhiBd);    
    
    % Half step y
    switch typeVel
        case 'att'
            Vy = fDirectionY(b,Vxo,Vyo,Ex,Ey,PhiDef); 
        otherwise
            Vy = Vyo;
    end
    btemp = fUpdateY(btemp,CFL,Vy, PhiDef, PhiAt, PhiBd);
    
    %% Total mass update
    % 03/07 - Correction of saving synchronisation
    B(n) = sum(sum(Dx*Dy*b));
   
    if ismember(n,nSaveVec)
        DtSave = Dt*(n)-TT(NSave+1);
        CFLtemp = DtSave/Dx;
%         Dt*(n-1)-DtSave
%         NSave+1

        % Half step x
        switch typeVel
            case 'att'
                Vx = fDirectionX(b,Vxo,Vyo,Ex,Ey,PhiDef); 
            otherwise
                Vx = Vxo;
        end
        b = fUpdateX(b,CFLtemp,Vx, PhiDef, PhiAt, PhiBd);

        % Half step y
        switch typeVel
            case 'att'
                Vy = fDirectionY(b,Vxo,Vyo,Ex,Ey,PhiDef); 
            otherwise
                Vy = Vyo;
        end
        b = fUpdateY(b,CFLtemp,Vy, PhiDef, PhiAt, PhiBd);
        
        % Save
        tsave =  Dt*(n-1)-DtSave;
        s = sprintf('%03s',num2str(NSave,'%d'));
        save([folderData '\' nameFile '-' s],...
            'X','Y','b','Dx','Dy','Vx','Vy','tsave')
        NSave = NSave + 1;
        
        % Estimation time
        plotTime(toc, n, Nt)
        pause(0.01)
    end
    b = btemp;
end


%% Post treatment
plotFinalTime();
Tcomp = toc;
tt = Dt*(1:Nt);
save([folderData '\' nameFile '-init'],...
    'B','Tcomp','tt','TT','-append')
% figure
% plotTotalMass(Dt*(1:Nt),B);

% savefig([nameFolder '\n' num2str(nx,'%d') '-' nameFile '.fig'])
% savefig([nameFig '\Fig_' s '.fig'])
% print([nameFig 'Png_' s],'-dpng')

