function [X,Y,Dx,Dy] = fGridGeneration(nx,Domain)
%fGridGeneration: Generate a grid [X,Y], definition regions for PhiBd,
%PhiAt for a Space in a Domain
% V1.1: Obstacle generation implemented
% V2.0: Downgrade to only grid generation and spacing

% Cells grid
Dx = (Domain(2)-Domain(1))/nx;
xx = linspace(Domain(1)+Dx/2,Domain(2)-Dx/2,nx);
Dy = Dx; ny = ceil((Domain(4)-Domain(3))/Dy);
yy = linspace(Domain(3)+Dy/2,Domain(4)-Dy/2,ny);
[X,Y] = meshgrid(xx,yy);

% % Definition regions
% PhiAt = ones(size(X))...
%     .*(X>Attractant(1)-Dx/2).*(X<Attractant(2)+Dx/2)...
%     .*(Y>Attractant(3)-Dy/2).*(Y<Attractant(4)+Dy/2);
% 
% PhiBd = ones(size(X)).*(X<Space(1)+Dx/2)...
%     +ones(size(X)).*(X>Space(2)-Dx/2)...
%     +ones(size(X)).*(Y<Space(3)+Dy/2).*(X<=Space(2)-Dx/2).*(X>=Space(1)+Dx/2)...
%     +ones(size(X)).*(Y>Space(4)-Dy/2).*(X<=Space(2)-Dx/2).*(X>=Space(1)+Dx/2);
% PhiBd = PhiBd.*(1-PhiAt);

% % Obstacle generation
% Quadrants
% s = size(Obstacle);
% No = s(1);
% for n = 1:No
%     PhiBd = PhiBd...
%     +ones(size(X)).*(X>=Obstacle(n,1)+Dx/2).*(X<=Obstacle(n,2)-Dx/2)...
%     .*(Y>=Obstacle(n,3)+Dy/2).*(Y<=Obstacle(n,4)-Dy/2);
% end

% % Discs
% s = size(Obstacle);
% No = s(1);
% for n = 1:No
% %     R = sqrt(Obstacle(n,2)-Obstacle(n,1)-Dx);
%     R = 0.5*(Obstacle(n,2)-Obstacle(n,1));
%     Xo = Obstacle(n,1)+R/2;
%     Yo = Obstacle(n,3)+R/2;
% %     [Xo-R^2, Xo+R^2, Yo-R^2, Yo+R^2]
%     if (Xo-2*sqrt(R))>Space(1) && (Xo+2*sqrt(R))<Space(2)...
%             && (Yo-2*sqrt(R))>Space(3)&& (Yo+2*sqrt(R))<Space(4)
% %     if 1
%         PhiBd = PhiBd...
%         +ones(size(X)).*(((X-Xo).^2+(Y-Yo).^2)<=R);
% %         contourf(X,Y,PhiBd)
% %         pause
%     end
% end
% 
% 
% % Limiter
% PhiBd = PhiBd.*(PhiBd<=1).*(PhiBd>=0)+1.*(PhiBd>1)+0.*(PhiBd<0);
% 
% PhiDef = ones(size(X))-PhiBd-PhiAt;
end

