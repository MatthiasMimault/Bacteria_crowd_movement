function [domAt, domBd, domDef, domSrc] = fRegionGeneration(...
X, Y, Space, Attractant, Source)
global Dx Dy typeAt typeObs typeSrc dev
% Define the numerical domain of definition, attractant, boundary and
% source based on geometric input and a grid X;
% Pending: obstacle generation

% Attractant region
domAt = fAttractantGen(X,Y,Dx,Dy,Attractant,typeAt);

% Source region
if dev 
    domSrc = fSourceGen_dev(X,Y,Dx,Dy,Source,Space,typeSrc);
else
    domSrc = fSourceGen(X,Y,Dx,Dy,Source,typeSrc);
end

% Boundary region
domBd = fBoundaryGen(X,Y,Dx,Dy,Space,Attractant,domAt,domSrc, typeObs);

% Definition domain
domDef = 1-domBd;
end

function domAt = fAttractantGen(X,Y,Dx,Dy,Attractant,typeAt)
% domAt = zeros(size(X));
switch typeAt
    case 'root'
        Xo = (Attractant(1)+Attractant(2))/2; Yo = Attractant(3);
        R = (Attractant(2)-Attractant(1))/2;
        domAt = ones(size(X))...
        .*(X>Attractant(1)-Dx/2).*(X<Attractant(2)+Dx/2)...
        .*(Y>Attractant(3)-Dy/2).*(Y<Attractant(4)+Dy/2)...
        +ones(size(X)).*(((X-Xo).^2+(Y-Yo).^2)<=R^2).*(Y<=Attractant(3)-Dy/2);
    case 'up'
    domAt = ones(size(X))...
        .*(X>Attractant(1)-Dx/2).*(X<Attractant(2)+Dx/2)...
        .*(Y>Attractant(3)-Dy/2).*(Y<Attractant(4)+Dy/2);
    case 'none'
        domAt = zeros(size(X));
    otherwise
    domAt = ones(size(X))...
        .*(X>Attractant(1)-Dx/2).*(X<Attractant(2)+Dx/2)...
        .*(Y>Attractant(3)-Dy/2).*(Y<Attractant(4)+Dy/2);
end
end

function domSrc = fSourceGen(X,Y,Dx,Dy,Source,typeSrc)
% domBd = zeros(size(X));
switch typeSrc
    case 'bottom'
        domSrc = (X>Source(1)-Dx/2).*(X<Source(2)+Dx/2)...
        .*(Y>Source(3)-Dy/2).*(Y<Source(4)+Dy/2);
    case 'Ulow'
        domSrc = (X>Source(1)-Dx/2).*(X<Source(2)+Dx/2)...
        .*(Y>Source(3)-Dy/2).*(Y<Source(4)+Dy/2);
        domSrc = domSrc.*(Y<0);    
    otherwise
        domSrc = zeros(size(X));
end
end

function domSrc = fSourceGen_dev(X,Y,Dx,Dy,Source,Space,typeSrc)
% domBd = zeros(size(X));
switch typeSrc
    case 'bottom'
        domSrc = (X>Source(1)-Dx/2).*(X<Source(2)+Dx/2)...
        .*(Y>Source(3)-Dy/2).*(Y<Source(4)+Dy/2);
    case 'Ulow'
        domSrc = (X>Source(1)-Dx/2).*(X<Source(2)+Dx/2)...
        .*(Y>Source(3)-Dy/2).*(Y<Source(4)+Dy/2);
        domDef = (X>Space(1)-Dx/2).*(X<Space(2)-Dx/2)...
        .*(Y>Space(3)-Dy/2).*(Y<Space(4)-Dy/2);
        domSrc = domSrc.*(1-domDef);    
    otherwise
        domSrc = zeros(size(X));
end
end

function domBd = fBoundaryGen(X,Y,Dx,Dy,Space,Attractant,domAt,domSrc, typeObs)
% Generate wall boundary, and obstacles. It includes domAt and domSrc
global Ro ro
%% Initialisation with wall boundary
% domBd = (X<Space(1)+Dx/2)+(X>Space(2)-Dx/2)...
%     +(Y<Space(3)+Dy/2).*(X<=Space(2)-Dx/2).*(X>=Space(1)+Dx/2)...
%     +(Y>Space(4)-Dy/2).*(X<=Space(2)-Dx/2).*(X>=Space(1)+Dx/2);
domBd = 1-(X>=Space(1)-Dx/2).*(X<=Space(2)-Dx/2).*(Y>=Space(3)-Dy/2).*(Y<=Space(4)-Dy/2);

%% Inclusion of Attractant and Source
domBd = max(domAt,domBd);
% Required: find a test of Attractant-Space-Source connectivity

%% Obstacle generation
% Required on V3.1
switch typeObs
    case 'particles'
        % Ro length between two consecutive centres, ro radius obstacle
%     Ro = 8; ro = 4;
    Obstacle = fObstacleGeneration(Space,Attractant,Ro,ro);
    s = size(Obstacle);
    No = s(1);
    
    for n = 1:No
    %     R = sqrt(Obstacle(n,2)-Obstacle(n,1)-Dx);
        R = 0.5*(Obstacle(n,2)-Obstacle(n,1));
        Xo = Obstacle(n,1)+Ro/2;
        Yo = Obstacle(n,3)+Ro/2;
        
        domBd = domBd+ones(size(X)).*(((X-Xo).^2+(Y-Yo).^2)<=R^2);
        
        
    end
    otherwise
        domBd = domBd.*(1-domSrc);
end
end
