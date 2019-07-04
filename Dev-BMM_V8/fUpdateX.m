function b = fUpdateX(b,CFL,Vx,PhiDef,PhiAt,PhiBd)
% CFL, PhiDef, PhiAt should be global

%% Parameters
s = size(b);
ny = s(2); % To be checked with non uniform grid

%% Flux
% % Hughes
% f = @(u) u.*(1-u);

% Linear
f = @(u) u;

%% Numerical scheme
% Lax Friedrichs
LxF = @(fp,fm,up,um) 0.5*(fp+fm-(up-um)/CFL);

%% Boundary condition
% Dirichlet 0
% Upx = [b,zeros(ny,1)];
% Umx = [zeros(ny,1),b];
% fpx = [-Vx,zeros(ny,1)].*f(Upx);
% fmx = [zeros(ny,1),-Vx].*f(Umx);

% Neumann 0 - Replicate
% PhiBdL2 = [zeros(ny,1),PhiBd(:,1:end-1)].*PhiDef;
% PhiBdR2 = [PhiBd(:,2:end),zeros(ny,1)].*PhiDef;
% UbdL = [b(:,2:end),zeros(ny,1)].*PhiBdL2
% UbdR = [zeros(ny,1),b(:,1:end-1)].*PhiBdR2
% Upx = [b+UbdR,zeros(ny,1)];
% Umx = [zeros(ny,1),b+UbdL];
% 
% Vpx = [Vx+[zeros(ny,1),Vx(:,1:end-1)].*PhiBdR2,zeros(ny,1)];
% Vmx = [zeros(ny,1),Vx+[Vx(:,2:end),zeros(ny,1)].*PhiBdL2];
% 
% fpx = -Vpx.*f(Upx)
% fmx = -Vmx.*f(Umx)

% Neumann 2 - Replicate U
PhiBdL = [PhiBd(:,1:end-1).*PhiDef(:,2:end),zeros(ny,1)];
PhiBdR = [zeros(ny,1),PhiBd(:,2:end).*PhiDef(:,1:end-1)];
UbdL = [b(:,2:end),zeros(ny,1)].*PhiBdL;
UbdR = [zeros(ny,1),b(:,1:end-1)].*PhiBdR;
Upx = [b+UbdR,zeros(ny,1)];
Umx = [zeros(ny,1),b+UbdL];

% Replicate and negate f
VbdL = [Vx(:,2:end),zeros(ny,1)].*PhiBdL;
VbdR = [zeros(ny,1),Vx(:,1:end-1)].*PhiBdR;
Vpx = [Vx-VbdR,zeros(ny,1)];
Vmx = [zeros(ny,1),Vx-VbdL];

fpx = Vpx.*f(Upx);
fmx = Vmx.*f(Umx);
    
%% Quantity update
LF = LxF(fpx,fmx,Upx,Umx);
b = b - CFL*(LF(:,2:end)-LF(:,1:end-1));

%% Boundary conditions 2
% Dirichlet 0 + flux limiter at exit
% b = b.*PhiDef+PhiAt.*min(b,0.2*ones(size(b)));
% Dirichlet 0
b = b.*PhiDef;
end

