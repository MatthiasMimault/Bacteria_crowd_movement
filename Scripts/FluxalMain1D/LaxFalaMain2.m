function [e,X,B0,BS,B1] = LaxFalaMain2(Dx,type)
%% Lax Friderichs a la main 2
% Normalised, tentative of convergence


% Initialisation and parameters
T = 1; X0 = 1;
a = 0.5; CFL = 0.9;
X = 0+Dx/2:Dx:4-Dx/2;
B0 = 0.5*(X>=X0);
BS = 0.5*(X>=X0+T*a);

% Time matter
Dt = CFL*Dx/a;
tt = 0:Dt:T;
nn = 1:length(tt)-1;

f = @(u) a*u;
switch type
    case 'Upwind'
        F = @(umm,um,up)  f(um);
    case 'Lax-Friedrichs'
        F = @(umm,um,up) 0.5*(f(um)+f(up)-Dx/Dt*(up-um));
    case 'Lax-Wendroff'
        F = @(umm,um,up) 0.5*(f(um)+f(up)-Dt/Dx*a*a*(up-um));
    case 'Lax-Wendroff2'
        F = @(umm,um,up)...
            (min(a,0)*up+max(a,0)*um)+0.5*abs(a)*(1-Dt/Dx*abs(a))*(up-um);        
    case 'Beam'
        F = @(umm,um,up) a*um+0.5*a*(1-Dt/Dx*a)*(um-umm);
    otherwise
        F = @(umm,um,up) 0.5*(f(um)+f(up)-Dx/Dt*(up-um));        
end

B = B0;

%% Loop
for n = nn
    BMM = [0, 0, B(1:end-2)];
    BM = [0 B(1:end-1)];
    BP = [B(2:end) 0.5];
    B = B-Dt/Dx*(F(BM,B,BP)-F(BMM,BM,B));
end

B1 = B;
e = Dx*sum(abs(BS-B1));

% figure
% plot(X,BS,X,B1)

% hold on
% stairs(X-Dx/2,BS)
% stairs(X-Dx/2,B1)
% stairs(X-Dx/2,B0)
% hold off