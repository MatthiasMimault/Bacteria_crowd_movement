function b = fDensityUpdate(Dt,b,Vx,Vy,Bdy,Src)
%FUPDATE updates the value of density b over a time step Dt with advection
%V, source Src and diffusion
global dev

typeDim = '2D';
% typeDim = 'none';

switch dev
    case 'stale'
        % Advection diffusion source
        % x direction
        b = fAdvectionX_stale(Dt,b,Vx,Bdy,Src);
        b = fDiffusionX_stale(Dt,b,Bdy,Src);

        % y direction
        b = fAdvectionY_stale(Dt,b,Vy,Bdy,Src);
        b = fDiffusionY_stale(Dt,b,Bdy,Src);
    case 'dev'
        switch typeDim
            case '1D'
                b = fAdvectionX_dev(Dt,b,Vx,Bdy,Src);
        %         b = fDiffusionX_dev(Dt,b,Bdy,Src);
            case '2D'
                b = fAdvectionX_dev(Dt,b,Vx,Bdy,Src);
                b = fAdvectionY_dev(Dt,b,Vy,Bdy,Src);
        end
    otherwise
        b = fAdvectionX_dev(Dt,b,Vx,Bdy,Src);
        b = fAdvectionY_dev(Dt,b,Vy,Bdy,Src);
end       
end

function b = fAdvectionX_stale(Dt,b,Vx,Bdy,Src)
% CFL, Def, PhiAt should be global
global Dx A 

%% Parameters
s = size(b);
ny = s(2); % To be checked with non uniform grid
% Def = 1-Bdy+Src;
Bdy = Bdy-Src;
Def = 1-Bdy-Src;
SrcValue = 0.5;

b = b+SrcValue*Src;

%% Flux
% % % Hughes
% f = @(u) A*u.*(1-u);

% Linear
f = @(u) A*u;

%% Numerical scheme
% Lax Friedrichs
LxF = @(fp,fm,up,um) 0.5*(fp+fm-Dx/Dt*(up-um));

%% Boundary condition
% Dirichlet 0
% Upx = [b,zeros(ny,1)];
% Umx = [zeros(ny,1),b];
% fpx = [-Vx,zeros(ny,1)].*f(Upx);
% fmx = [zeros(ny,1),-Vx].*f(Umx);

% Neumann 0 - Replicate
BdyL = [Bdy(:,1:end-1).*Def(:,2:end),zeros(ny,1)];
BdyR = [zeros(ny,1),Bdy(:,2:end).*Def(:,1:end-1)];
UbdL = [b(:,2:end),zeros(ny,1)].*BdyL;
UbdR = [zeros(ny,1),b(:,1:end-1)].*BdyR;
Upx = [b+UbdR,zeros(ny,1)];
Umx = [zeros(ny,1),b+UbdL];

% Replicate and negate f
VbdL = [Vx(:,2:end),zeros(ny,1)].*BdyL;
VbdR = [zeros(ny,1),Vx(:,1:end-1)].*BdyR;
Vpx = [Vx-VbdR,zeros(ny,1)];
Vmx = [zeros(ny,1),Vx-VbdL];

fpx = Vpx.*f(Upx);
fmx = Vmx.*f(Umx);
    
%% Quantity update
LF = LxF(fpx,fmx,Upx,Umx);
b = b - Dt/Dx*(LF(:,2:end)-LF(:,1:end-1));

%% Boundary conditions 2
% Dirichlet 0 + flux limiter at exit
% b = b.*Def+PhiAt.*min(b,0.2*ones(size(b)));
% Dirichlet 0
b = b.*Def;
end

function b = fAdvectionY_stale(Dt,b,Vy,Bdy,Src)
global A Dy BactValue

%% Parameters
s = size(b);
nx = s(1); % To be checked with non uniform grid
% Def = 1-Bdy+Src;
Bdy = Bdy-Src;
Def = 1-Bdy-Src;
% SrcValue = 0.5;

b = b+BactValue*Src;

%% Flux
% % Hughes
% f = @(u) A*u.*(1-u);

% Linear
f = @(u) A*u;

%% Numerical scheme
% Lax Friedrichs
LxG = @(fp,fm,up,um) 0.5*(fp+fm-Dy/Dt*(up-um));

%% Boundary condition
% Dirichlet 0
% Upy = [b;zeros(1,nx)];
% Umy = [zeros(1,nx);b];
% fpy = [-Vy;zeros(1,nx)].*f(Upy);
% fmy = [zeros(1,nx);-Vy].*f(Umy);

% % Neumann 0 - Replicate
BdyD = [Bdy(1:end-1,:).*Def(2:end,:);zeros(1,nx)];
BdyU = [zeros(1,nx);Bdy(2:end,:).*Def(1:end-1,:)];
UbdD = [b(2:end,:);zeros(1,nx)].*BdyD;
UbdU = [zeros(1,nx);b(1:end-1,:)].*BdyU;
Upy = [b+UbdU;zeros(1,nx)];
Umy = [zeros(1,nx);b+UbdD];

% Replicate and negate f
VbdD = [Vy(2:end,:);zeros(1,nx)].*BdyD;
VbdU = [zeros(1,nx);Vy(1:end-1,:)].*BdyU;
Vpy = [Vy-VbdU;zeros(1,nx)];
Vmy = [zeros(1,nx);Vy-VbdD];

fpy = Vpy.*f(Upy);
fmy = Vmy.*f(Umy);
% figure
% contourf(fpy)
% figure
% contourf(fmy)
% pause
    
%% Quantity update
LG = LxG(fpy,fmy,Upy,Umy);
b = b - Dt/Dy*(LG(2:end,:)-LG(1:end-1,:));

%% Boundary conditions 2
% Dirichlet 0 + flux limiter at exit
% b = b.*Def+PhiAt.*min(b,0.2*ones(size(b)));
% Dirichlet 0
b = b.*Def;
end

function b = fDiffusionX_stale(Dt,b,Bdy,Src)
global Dx C BactValue

%% Parameters
s = size(b);
ny = s(2);
Bdy = Bdy-Src;
Def = 1-Bdy-Src;
% SrcValue = 0.5;

b = b+BactValue*Src;

%% Neumann boundary
bp = [zeros(ny,1) b(:,1:end-1)].*Bdy+b;
bm = [b(:,2:end) zeros(ny,1)].*Bdy+b;

b = b+Dt*C/Dx/Dx*(bm-2*b+bp);
b = b.*Def;
end

function b = fDiffusionY_stale(Dt,b,Bdy,Src)
%% Parameters
global Dy C
s = size(b);
nx = s(1);
Bdy = Bdy-Src;
Def = 1-Bdy-Src;
SrcValue = 0.5;

b = b+SrcValue*Src;

%% Neumann boundary
bu = [zeros(1,nx);b(1:end-1,:)].*Bdy+b;
bd = [b(2:end,:);zeros(1,nx)].*Bdy+b;

b = b+Dt/Dy/Dy*C*(bd-2*b+bu);
b = b.*Def;
end

function b = fAdvectionX_dev(Dt,b,Vx,Bdy,Src)
% CFL, Def, PhiAt should be global
global Dx A D BactValue

%% Parameters
s = size(b);
ny = s(1); % To be checked with non uniform grid
Bdy = Bdy-Src;
Def = 1-Bdy-Src;
% SrcValue = 0.5;

b = b+BactValue*Src;

%% Flux
% % % Hughes
% f = @(u) A*u.*(1-u);

% Linear
f = @(u) A*u;

%% Numerical scheme
% Lax Friedrichs
LxF =  @(fp,fm,up,um) 0.5*(fp+fm);
% LxFr = @(fp,fm,up,um) 0.5*(fp+fm-A*(up-um));

%% Boundary condition
% Neumann 0 - Replicate
BdyL = [Bdy(:,1:end-1).*Def(:,2:end),zeros(ny,1)];
BdyR = [zeros(ny,1),Bdy(:,2:end).*Def(:,1:end-1)];
UbdL = [b(:,2:end),zeros(ny,1)].*BdyL;
UbdR = [zeros(ny,1),b(:,1:end-1)].*BdyR;
Upx = [b+UbdR,zeros(ny,1)];
Umx = [zeros(ny,1),b+UbdL];

% Replicate and negate f
VbdL = [Vx(:,2:end),zeros(ny,1)].*BdyL;
VbdR = [zeros(ny,1),Vx(:,1:end-1)].*BdyR;
Vpx = [Vx-VbdR,zeros(ny,1)];
Vmx = [zeros(ny,1),Vx-VbdL];

fpx = Vpx.*f(Upx);
fmx = Vmx.*f(Umx);
    
%% Quantity update
LF = LxF(fpx,fmx,Upx,Umx);
DF = Upx-Umx;

b = b - Dt/Dx*(LF(:,2:end)-LF(:,1:end-1))...
    +0.5*Dt*max(A/Dx,2*D/Dx/Dx)*(DF(:,2:end)-DF(:,1:end-1));

%% Boundary conditions 2
% Dirichlet 0 + flux limiter at exit
% b = b.*Def+PhiAt.*min(b,0.2*ones(size(b)));
% Dirichlet 0
b = b.*Def;
end

function b = fAdvectionY_dev(Dt,b,Vy,Bdy,Src)
% CFL, Def, PhiAt should be global
global Dx A D BactValue

%% Parameters
s = size(b);
nx = s(2); % 
Bdy = Bdy-Src;
Def = 1-Bdy-Src;
b = b+BactValue*Src;

%% Flux
% % Hughes
% f = @(u) A*u.*(1-u);

% Linear
f = @(u) A*u;

%% Numerical scheme
% Lax Friedrichs
LxG = @(fp,fm,up,um) 0.5*(fp+fm);

%% Boundary condition
% % Neumann 0 - Replicate
BdyD = [Bdy(1:end-1,:).*Def(2:end,:);zeros(1,nx)];
BdyU = [zeros(1,nx);Bdy(2:end,:).*Def(1:end-1,:)];
UbdD = [b(2:end,:);zeros(1,nx)].*BdyD;
UbdU = [zeros(1,nx);b(1:end-1,:)].*BdyU;
Upy = [b+UbdU;zeros(1,nx)];
Umy = [zeros(1,nx);b+UbdD];

% Replicate and negate f
VbdD = [Vy(2:end,:);zeros(1,nx)].*BdyD;
VbdU = [zeros(1,nx);Vy(1:end-1,:)].*BdyU;
Vpy = [Vy-VbdU;zeros(1,nx)];
Vmy = [zeros(1,nx);Vy-VbdD];

fpy = Vpy.*f(Upy);
fmy = Vmy.*f(Umy);
    
%% Quantity update
LG = LxG(fpy,fmy,Upy,Umy);
DG = Upy-Umy;

b = b - Dt/Dx*(LG(2:end,:)-LG(1:end-1,:))...
    +0.5*Dt*max(A/Dx,2*D/Dx/Dx)*(DG(2:end,:)-DG(1:end-1,:));

%% Boundary conditions 2
% Dirichlet 0
b = b.*Def;
end




