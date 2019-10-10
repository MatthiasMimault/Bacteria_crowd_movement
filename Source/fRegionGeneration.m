function [domAt, domBd, domDef, domSrc] = fRegionGeneration(...
X, Y, Domain, Space, Attractant, Source)
global Dx Dy typeAt typeSrc
% Define the numerical domain of definition, attractant, boundary and
% source based on geometric input and a grid X;

% Attractant region
domAt = fAttractantGen(X,Y,Dx,Dy,Attractant,typeAt);

% Source region
domSrc = fSourceGen(X,Y,Dx,Dy,Source,typeSrc);

% Boundary region
domBd = zeros(size(X));

domDef = zeros(size(X));
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
        domSrc = ones(size(X))...
        -(X>Source(1)-Dx/2).*(X<Source(2)+Dx/2)...
        .*(Y>Source(3)-Dy/2).*(Y<Source(4)+Dy/2);
        domSrc = domSrc.*(Y<0);
    
    otherwise
        domSrc = zeros(size(X));
end
end

function domBd = fBoundaryGen(X,Y,Dx,Dy,Space,domAt,domSrc)

end