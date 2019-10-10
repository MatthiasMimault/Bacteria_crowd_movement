function b = fUpdateY(b,Dt,Dy,A,Vy,PhiDef,PhiAt,PhiBd)
% CFL, PhiDef, PhiAt should be global

%% Parameters
s = size(b);
nx = s(1); % To be checked with non uniform grid

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
PhiBdD = [PhiBd(1:end-1,:).*PhiDef(2:end,:);zeros(1,nx)];
PhiBdU = [zeros(1,nx);PhiBd(2:end,:).*PhiDef(1:end-1,:)];
UbdD = [b(2:end,:);zeros(1,nx)].*PhiBdD;
UbdU = [zeros(1,nx);b(1:end-1,:)].*PhiBdU;
Upy = [b+UbdU;zeros(1,nx)];
Umy = [zeros(1,nx);b+UbdD];

% Replicate and negate f
VbdD = [Vy(2:end,:);zeros(1,nx)].*PhiBdD;
VbdU = [zeros(1,nx);Vy(1:end-1,:)].*PhiBdU;
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
% b = b.*PhiDef+PhiAt.*min(b,0.2*ones(size(b)));
% Dirichlet 0
% b = b.*PhiDef;
end

