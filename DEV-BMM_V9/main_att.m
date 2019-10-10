%% Main file for Bacterial movement model
% V1.70 - Attractive root, repulsive region, weno3
% V1.71 - Tentative automatic file naming aborted, correction of NsaveVect
% generation (validated mostly, no wrong case found for CFL 0.1-0.8)
% V1.72 - Correction of behaviour at boundary: Solved, but Velocity
% orientation troubling
% V1.8 - Stable method. Problem at vector field generation
close all
clear

%% Parameters
% Settings: Time T, Advection velocity A, Attractant trail B
% Obstacles: ro is the size (side) of a square obstacle and
% Ro the spacing between the center of two obstacles
% R is the radius of interaction for the perceived density
T = 10; A = 0.03; B = 1; 
R = 0.5;
Nfiles = 10; 
typeTest = 'CLB-Rdim';
dateTest = '0907';
runnb = '1';
% A standstill, B advection, C adv-obstacles, D attraction, E att-obstacles
% F Root alone
type = 'D'; 

% Numerical parameters
nx = 1200;
CFL = 0.9;
ny = nx;

% nameFile = fCaseNaming(type,suffix,nx,nameFolder);

% Domain
Domain = [-25, 25, -25, 25];
Space = [-20, 20, -20, 20];
Attractant = [-1 1 10 20];
% Attractant = [-20 20 20 25];

% Bacteria
InitSpace = [-19, 19, -19, -5];
InitValue = 0.1;

switch type
    case 'A'
        suffix = 'Standstill';
        typeAt = 'none';
        typeObs = 'none';
        typeVel = 'stand';
    case 'B'
        suffix = 'Advection';
        typeAt = 'up';
        typeObs = 'none';
        typeVel = 'adv';
    case 'C'
        suffix = 'AdvObstacles';
        typeAt = 'up';
        typeObs = 'particles';
        typeVel = 'adv';
    case 'D'
        suffix = 'Attraction';
        typeAt = 'up';
        typeObs = 'none';
        typeVel = 'att';
    case 'E'
        suffix = 'AttObstacles';
        typeAt = 'up';
        typeObs = 'particles';
        typeVel = 'att';
    case 'F'
        suffix = 'AdvRoot';
        typeAt = 'root';
        typeObs = 'none';
        typeVel = 'adv';
    case 'Fb'
        suffix = 'StdRoot';
        typeAt = 'none';
        typeObs = 'none';
        typeVel = 'stand';
    case 'G'
        suffix = 'ObsRoot';
        typeAt = 'root';
        typeObs = 'particles';
        typeVel = 'adv';
    case 'H'
        suffix = 'ObsAttRoot';
        typeAt = 'root';
        typeObs = 'particles';
        typeVel = 'att';
        
    otherwise
        disp('Unknown case. Abort')
        return
end
nameFolder = [typeTest '-' dateTest];
nameFile = [type runnb '-' suffix];


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
folderData = [nameFolder '\Data-' nameFile '-' num2str(nx)];
mkdir(nameFolder)
mkdir(folderData)
addpath('..\Toolbox')  

% Define regions
% Domain generation (2D) (Finite Volumes Cartesian grid)
% V2.0 Generation ok, to be tested. To move Dxy inside function
% V3.0 Breakdown of region creation with string options
% [X,Y,Dx,Dy,PhiBd,PhiAt,PhiDef] = ...
%     fGridGeneration(nx,ny, Domain, Space, Attractant, Obstacles);
% V1.81 - Adaptation to root space which at the same time Bd and At

[X,Y,Dx,Dy] = fGridGeneration(nx,ny,Domain);
PhiAt = fAttractantGeneration(X,Y,Attractant,typeAt);
PhiBd = fBoundaryGeneration(X,Y,PhiAt,Space,Attractant,typeAt,typeObs);
PhiDef = ones(size(X))-min(1,PhiAt+PhiBd);

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
save([folderData '\' nameFile '-' num2str(nx) '-init'],'PhiDef','PhiBd','PhiAt',...
    'Nt','dN','T','Nfiles', 'Axis', 'nameFolder','nameFile','nx')
save([folderData '\' nameFile '-' num2str(nx) '-000'])

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
%     plotSurf(X,Y,b,0,0,Axis)
%     max(max(b))
%     min(min(b))
%     pause
   
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
        save([folderData '\' nameFile '-' num2str(nx) '-' s],...
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
save([folderData '\' nameFile '-' num2str(nx) '-init'],...
    'B','Tcomp','tt','TT','-append')
% figure
% plotTotalMass(Dt*(1:Nt),B);

% savefig([nameFolder '\n' num2str(nx,'%d') '-' nameFile '.fig'])
% savefig([nameFig '\Fig_' s '.fig'])
% print([nameFig 'Png_' s],'-dpng')

