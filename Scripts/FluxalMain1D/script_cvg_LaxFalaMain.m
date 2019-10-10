%% Convergence script for the function Lax Friedrichs a la main

% Vector Dx
N = 12;
DX = 0.125./(2.^(1:N));

% Vector e
E = zeros(1,length(DX));

%% Loop
for n=1:N
    if (n>=11)
        n
    end
    [E(n),X,B0,BS,B1] = AdvectionAbs(DX(n),'Lax-Wendroff');
%     [DX(n),E(n)]
end

% cvg_order = log(abs(E(:,1:end-1)-E(:,2:end)))./log(abs(DX(1:end-1)-DX(2:end)))
cvg_order = log2(E(:,1:end-1)./E(:,2:end))
loglog(DX,E)
plot(X,B1,X,BS)