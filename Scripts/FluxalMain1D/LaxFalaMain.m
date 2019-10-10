function [e,X,B0,BS,B1] = LaxFalaMain(Dx)
%% Lax Friderichs a la main

% Initialisation and parameters
T = 15;
a = 0.03; CFL = 0.9;
X = 0+Dx/2:Dx:2-Dx/2;
B0 = 0.5*(X>=1);
BS = 0.5*(X>=1+T*a);

% Time matter
Dt = CFL*Dx/a;
tt = 0:Dt:T;
nn = 1:length(tt)-1;



f = @(u) a*u;
F = @(um,up) 0.5*(f(um)+f(up)-Dx/Dt*(up-um));

B = B0;

%% Loop
for n = nn
    BM = [0 B(1:end-1)];
    BP = [B(2:end) 0.5];
    B = B-Dt/Dx*(F(B,BP)-F(BM,B));
end

B1 = B;
e = Dx*sum(abs(BS-B1));

hold on
% stairs(X-Dx/2,BS)
% stairs(X-Dx/2,B1)
stairs(X-Dx/2,B0)
hold off