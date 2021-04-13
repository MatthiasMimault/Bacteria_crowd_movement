function b = fDensityUpdate(Dt,b,Vx,Vy,Bdy,Src)
%FUPDATE updates the value of density b over a time step Dt with advection
%V, source Src and diffusion
global dev

% typeDim = '2D';
% typeDim = 'none';

switch dev
    otherwise
        b = fAdvectionX(Dt,b,Vx);
        b = fAdvectionY(Dt,b,Vy);
end       
end

function b = fAdvectionX(Dt,b,Vx)
% CFL, Def, PhiAt should be global
global Dx A D BactValue domBd domSrc domAt domDef rateCol

%% Parameters
s = size(b);
ny = s(1);
Bdy = domBd-domAt;

b = b.*(1-domSrc)+BactValue*domSrc;

%% Flux
% % % Hughes
% f = @(u) A*u.*(1-u);

% Linear
f = @(u) A*u;
fm = @(u) -A*u;
fa = @(u) rateCol*u;

%% Numerical scheme
% Lax Friedrichs
LxF =  @(fp,fm,up,um) 0.5*(fp+fm);

%% Boundary condition
% Neumann 0 - Replicate
BdyL = [Bdy(:,1:end-1).*domDef(:,2:end),zeros(ny,1)];
BdyR = [zeros(ny,1),Bdy(:,2:end).*domDef(:,1:end-1)];
UbdL = [b(:,2:end),zeros(ny,1)].*BdyL;
UbdR = [zeros(ny,1),b(:,1:end-1)].*BdyR;

% Neumann surface At
AtL = [domAt(:,1:end-1).*domDef(:,2:end),zeros(ny,1)];
AtR = [zeros(ny,1),domAt(:,2:end).*domDef(:,1:end-1)];
UaL = [b(:,2:end),zeros(ny,1)].*AtL;
UaR = [zeros(ny,1),b(:,1:end-1)].*AtR;
Upx = [b+UbdR+UaR,zeros(ny,1)];
Umx = [zeros(ny,1),b+UbdL+UaL];

% Replicate and negate V
VbdL = [Vx(:,2:end),zeros(ny,1)].*BdyL;
VbdR = [zeros(ny,1),Vx(:,1:end-1)].*BdyR;
VaL = [Vx(:,2:end),zeros(ny,1)].*AtL;
VaR = [zeros(ny,1),Vx(:,1:end-1)].*AtR;
Vpx = [Vx+VbdR+VaR,zeros(ny,1)];
Vmx = [zeros(ny,1),Vx+VbdL+VaL];

% Flux and colonisation flux
% Fpx = f(Upx)+fo(Uapx);
% Fmx = f(Umx)+fo(Uamx);

fpx = Vpx.*(f(Upx).*[zeros(ny,1),domDef]...
    +fm(Upx).*[zeros(ny,1),domBd]...
    +fa(Upx).*[zeros(ny,1),domAt]);
fmx = Vmx.*(f(Umx).*[domDef,zeros(ny,1)]...
    +fm(Umx).*[domBd,zeros(ny,1)]...
    +fa(Umx).*[domAt,zeros(ny,1)]);
    
%% Quantity update
LF = LxF(fpx,fmx,Upx,Umx);
DF = Upx-Umx;

b = b - Dt/Dx*(LF(:,2:end)-LF(:,1:end-1))...
    +0.5*Dt*max(A/Dx,2*D/Dx/Dx)*(DF(:,2:end)-DF(:,1:end-1));

%% Boundary conditions 2
% Dirichlet 0
b = b.*(domDef-domSrc);
end

function b = fAdvectionY(Dt,b,Vy)
% CFL, Def, PhiAt should be global
global Dx A D BactValue domBd domSrc domAt domDef rateCol

%% Parameters
s = size(b);
nx = s(2); % 
Bdy = domBd-domAt;

b = b.*(1-domSrc)+BactValue*domSrc;

%% Flux
% % % Hughes
% f = @(u) A*u.*(1-u);

% Linear
f = @(u) A*u;
fm = @(u) -A*u;
fa = @(u) rateCol*u;

%% Numerical scheme
% Lax Friedrichs
LxG = @(fp,fm,up,um) 0.5*(fp+fm);

%% Boundary condition
% % Neumann 0 - Replicate
BdyD = [Bdy(1:end-1,:).*domDef(2:end,:);zeros(1,nx)];
BdyU = [zeros(1,nx);Bdy(2:end,:).*domDef(1:end-1,:)];
UbdD = [b(2:end,:);zeros(1,nx)].*BdyD;
UbdU = [zeros(1,nx);b(1:end-1,:)].*BdyU;

% Neumann surface At
AtD = [domAt(1:end-1,:).*domDef(2:end,:);zeros(1,nx)];
AtU = [zeros(1,nx);domAt(2:end,:).*domDef(1:end-1,:)];
UaD = [b(2:end,:);zeros(1,nx)].*AtD;
UaU = [zeros(1,nx);b(1:end-1,:)].*AtU;

Upy = [b+UbdU+UaU;zeros(1,nx)];
Umy = [zeros(1,nx);b+UbdD+UaD];

% Replicate and negate f
VbdD = [Vy(2:end,:);zeros(1,nx)].*BdyD;
VbdU = [zeros(1,nx);Vy(1:end-1,:)].*BdyU;
VaD = [Vy(2:end,:);zeros(1,nx)].*AtD;
VaU = [zeros(1,nx);Vy(1:end-1,:)].*AtU;
Vpy = [Vy-VbdU+VaU;zeros(1,nx)];
Vmy = [zeros(1,nx);Vy-VbdD+VaD];

% Flux and colonisation flux
fpy = Vpy.*(f(Upy).*[zeros(1,nx);domDef]...
    +fm(Upy).*[zeros(1,nx);domBd]...
    +fa(Upy).*[zeros(1,nx);domAt]);
fmy = Vmy.*(f(Umy).*[domDef;zeros(1,nx)]...
    +fm(Umy).*[domBd;zeros(1,nx)]...
    +fa(Umy).*[domAt;zeros(1,nx)]);
    
%% Quantity update
LG = LxG(fpy,fmy,Upy,Umy);
DG = Upy-Umy;

b = b - Dt/Dx*(LG(2:end,:)-LG(1:end-1,:))...
    +0.5*Dt*max(A/Dx,2*D/Dx/Dx)*(DG(2:end,:)-DG(1:end-1,:));

%% Boundary conditions 2
% Dirichlet 0
b = b.*(domDef-domSrc);
end
