%% Convergence order script
% V2.0 - Time evolution of convergence
% V2.1 - Compared to initial data (for standstill)

% nx = [150, 300, 600, 1200, 2400];
nx = [100, 200, 400, 800, 1600];
nf = 4;

nTests = length(nx);
e = zeros(nf, length(nx));
ddx = zeros(1, length(nx));

for n = 1:nTests
    folderData = [nameFolder '\Data-' nameFile '-' num2str(nx(n),'%d')];
    load([folderData '\' nameFile '-' num2str(nx(n),'%d') '-init'], 'Nt',...
        'dN', 'Axis');
    % load b0
    s = sprintf('000');
    load([folderData '\' nameFile '-' num2str(nx(n),'%d') '-' s],'b','Dx');
    b0 = b;
    ddx(n) = Dx;
    for m = 1:nf 
        
        % Load b2
        s = sprintf('%03s',num2str(m,'%d'));
        load([folderData '\' nameFile '-' num2str(nx(n),'%d') '-' s],...
            'b','Dx','Dy');
        
        % Error computation
        e(m,n) = Dx*Dy*sum(sum(abs(b0-b)));
%         contourf(abs(b0-b))
%         colorbar
%         pause
        pause(0.01)
    end
    
%     [e2, dx2] = plotErrorEstimate(nameFolder,nameFile,nx(n+1),nx(n+2));
%     p(n) = log(e2/e1)/log(dx2/dx1);
end

% Convergence rate, gamma, p, order of accuracy
% cvg_rate = log(e(:,1:end-1)./e(:,2:end))

e
% Convergence order, convergence speed for discretization methods
% cvg_order = log(abs(e(:,1:end-1)-e(:,2:end)))./log(abs(ddx(1:end-1)-ddx(2:end)))
cvg_order = log2(e(:,1:end-1)./e(:,2:end))

plot(1:nf,cvg_order)
legend('100','200','400','800')
title(['Cvg - ' nameFolder])
