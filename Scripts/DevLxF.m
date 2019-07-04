%% Study of Lax Friedrich scheme

%% Parameters
CFL = 0.5;
Dx = 0.001;
T = 1;
a = 1;
u0 = [0,0.5];

%% Variables
xx = -1+Dx/2:Dx:1-Dx/2;
uu = u0(1)*(xx<0)+u0(2)*(xx>0);
Dt = CFL*Dx;
tt = 0:Dt:T;

%% Functions
% f = @(u) a*u;
f = @(u) u.*(1-u);
% V0.0 LxF = @(fp,fm,up,um) 0.5*(fp+fm-0.5*Dx/Dt*(up-um));
LxF = @(fp,fm,up,um) 0.5*(fp+fm-Dx/Dt*(up-um));

%% Loop
for t=tt
    ff = LxF(f([uu,0]),f([0, uu]),[uu,0],[0, uu]);
    uu = uu - Dt/Dx*(ff(2:end)-ff(1:end-1));
    plot(xx,uu)
    axis([-1,1,0,1])
    pause(0.01)
end