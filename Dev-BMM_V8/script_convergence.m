%% Convergence order script
% V2.0 - Time evolution of convergence

% nx = [100, 200, 400, 800, 1600, 3200];
% nx = [100, 200, 400, 800, 1600];
% nx = [150, 300, 600, 1200, 2400];
nx = [150, 300, 600, 1200, 2400];
nf = 10;

nTests = length(nx);
e = zeros(nf, length(nx)-1);
ddx = zeros(1, length(nx)-1);

% for n = 1:nTests-2
%     [e1, dx1] = plotErrorEstimate(nameFolder,nameFile,nx(n),nx(n+1));
%     [e2, dx2] = plotErrorEstimate(nameFolder,nameFile,nx(n+1),nx(n+2));
%     p(n) = log(e2/e1)/log(dx2/dx1);
% end
% 
% nameFileTemp = nameFile;
% 
% for  n = nTests:-1:1
%     nameFileTemp = erase(nameFileTemp, ['-' num2str(nx(n))])
% end
% 
% nameFileTemp

for n = 1:nTests-1
    [e(:,n), ddx(n)] = plotErrorEstimate(nameFolder,nameFile,nx(n),nx(n+1));
    
%     [e2, dx2] = plotErrorEstimate(nameFolder,nameFile,nx(n+1),nx(n+2));
%     p(n) = log(e2/e1)/log(dx2/dx1);
end

% Convergence rate, gamma, p, order of accuracy
% cvg_rate = log(e(:,1:end-1)./e(:,2:end))


% Convergence order, convergence speed for discretization methods
cvg_order = log(abs(e(:,1:end-1)-e(:,2:end)))./log(abs(ddx(1:end-1)-ddx(2:end)))

plot(cvg_order)
legend('150','300','600')
title('Cvg - Advection')
