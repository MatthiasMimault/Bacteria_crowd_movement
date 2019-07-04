function b = fUpdateY(b,CFL,Vy,PhiDef,PhiAt,PhiBd)
% CFL, PhiDef, PhiAt should be global

%% Parameters
s = size(b);
nx = s(1); % To be checked with non uniform grid

%% Flux
% % Hughes
% f = @(u) u.*(1-u);

% Linear
f = @(u) u;

%% Numerical scheme
% Lax Friedrichs
LxG = @(fp,fm,up,um) 0.5*(fp+fm-1/CFL*(up-um));

%% Boundary condition
% Dirichlet 0
% Upy = [b;zeros(1,nx)];
% Umy = [zeros(1,nx);b];
% fpy = [-Vy;zeros(1,nx)].*f(Upy);
% fmy = [zeros(1,nx);-Vy].*f(Umy);

% % Neumann 0 - Replicate
PhiBdD = [PhiBd(1:end-1,:).*PhiDef(2:end,:);zeros(1,nx)];
PhiBdU = [zeros(1,nx);PhiBd(2:end,:).*PhiDef(1:end-1,:)];
UbdD = [b(2:end,:);zeros(1,nx)].*PhiBdD;
UbdU = [zeros(1,nx);b(1:end-1,:)].*PhiBdU;
Upy = [b+UbdU;zeros(1,nx)];
Umy = [zeros(1,nx);b+UbdD];

VbdD = [Vy(2:end,:);zeros(1,nx)].*PhiBdD;
VbdU = [zeros(1,nx);Vy(1:end-1,:)].*PhiBdU;
Vpy = [Vy-VbdU;zeros(1,nx)];
Vmy = [zeros(1,nx);Vy-VbdD];

fpy = Vpy.*f(Upy);
fmy = Vmy.*f(Umy);
% 
% % Neumann 2 - Replicate U
% PhiBdL = [PhiBd(:,1:end-1).*PhiDef(:,2:end),zeros(ny,1)];
% PhiBdR = [zeros(ny,1),PhiBd(:,2:end).*PhiDef(:,1:end-1)];
% UbdL = [b(:,2:end),zeros(ny,1)].*PhiBdL;
% UbdR = [zeros(ny,1),b(:,1:end-1)].*PhiBdR;
% Upx = [b+UbdR,zeros(ny,1)];
% Umx = [zeros(ny,1),b+UbdL];
% 
% % Replicate and negate f
% VbdL = [Vx(:,2:end),zeros(ny,1)].*PhiBdL;
% VbdR = [zeros(ny,1),Vx(:,1:end-1)].*PhiBdR;
% Vpx = [Vx-VbdR,zeros(ny,1)];
% Vmx = [zeros(ny,1),Vx-VbdL];
% 
% fpx = Vpx.*f(Upx);
% fmx = Vmx.*f(Umx);
    
%% Quantity update
LG = LxG(fpy,fmy,Upy,Umy);
b = b - CFL*(LG(2:end,:)-LG(1:end-1,:));

%% Boundary conditions 2
% Dirichlet 0 + flux limiter at exit
% b = b.*PhiDef+PhiAt.*min(b,0.2*ones(size(b)));
% Dirichlet 0
b = b.*PhiDef;
end

