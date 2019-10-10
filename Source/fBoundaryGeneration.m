function PhiBd = fBoundaryGeneration(X,Y,PhiAt,Space,Attractant,typeAt,typeObs)

%% Initialisation
Dx = X(1,2)-X(1,1);
Dy = Y(2,1)-Y(1,1);

%% Generation domain boundary
if strcmp(typeAt,'root')
%     Space = [-1, 1, -1, 1];
    PhiBd = ones(size(X))...
        -(X>Space(1)-Dx/2).*(X<Space(2)+Dx/2)...
        .*(Y>Space(3)-Dy/2).*(Y<Space(4)+Dy/2);
    PhiBd = PhiBd | PhiAt;
else
%     Space = [-1, 1, -1, 1];
    PhiBd = ones(size(X)).*(X<Space(1)+Dx/2)...
    +ones(size(X)).*(X>Space(2)-Dx/2)...
    +ones(size(X)).*(Y<Space(3)+Dy/2).*(X<=Space(2)-Dx/2).*(X>=Space(1)+Dx/2)...
    +ones(size(X)).*(Y>Space(4)-Dy/2).*(X<=Space(2)-Dx/2).*(X>=Space(1)+Dx/2);
    PhiBd = PhiBd.*(1-PhiAt);
end

%% Generation obstacles
if strcmp(typeObs,'particles')
    Ro = 10; ro = 4;
    Obstacle = fObstacleGeneration(Space,Attractant,Ro,ro);
    s = size(Obstacle);
    No = s(1);
    for n = 1:No
    %     R = sqrt(Obstacle(n,2)-Obstacle(n,1)-Dx);
        R = 0.5*(Obstacle(n,2)-Obstacle(n,1));
        Xo = Obstacle(n,1)+R/2;
        Yo = Obstacle(n,3)+R/2;
    %     [Xo-R^2, Xo+R^2, Yo-R^2, Yo+R^2]
        if ((Xo-2*sqrt(R))>Space(1) && (Xo+2*sqrt(R))<Space(2)...
                && (Yo-2*sqrt(R))>Space(3)&& (Yo+2*sqrt(R))<Space(4))
    %     if 1
            PhiBd = PhiBd...
            +ones(size(X)).*(((X-Xo).^2+(Y-Yo).^2)<=R^2);
    %         contourf(X,Y,PhiBd)
    %         pause
        end
    end
end
%     contourf(PhiBd)
%     pause

%     PhiBd = PhiBd.*(1-PhiAt);