function [Vxo, Vyo] = fVelocityGeneration(X, Y, domAt, domBd, domSrc)
%FVELOCITYGENERATION  generates a quiver vector field based on the domains
%of definition, the grid and type
global Dx Dy typeVel

domDef = 1-domBd+domSrc;

switch typeVel
    case 'adv-up'
        Vxo = zeros(size(X));
        Vyo = ones(size(X));
    case 'adv-right'
        Vxo = ones(size(X));
        Vyo = zeros(size(X));
    case {'adv-src', 'att-src'}
        Distance = fEikonalCost(X,Y,ones(size(X)),domBd-domAt-domSrc,domAt);
        [Vxo,Vyo] = fDiffFlex(Distance,domDef,Dx,Dy);
        Vn = sqrt(Vxo.^2+Vyo.^2+1-domDef);
        Vxo = -Vxo./Vn;
        Vxo(isnan(Vxo)) = 0;
        Vxo = Vxo.*(domDef);
        Vyo = -Vyo./Vn;
        Vyo(isnan(Vyo)) = 0;
        Vyo = Vyo.*(domDef);
    case {'att','att-adv'}
        PhiC = fEikonalCost(X,Y,ones(size(X)),domBd-domAt,domAt);
        [Vxo,Vyo] = fDiffFlex(PhiC,domDef,Dx,Dy);
        Vn = sqrt(Vxo.^2+Vyo.^2+domBd+domAt);
        Vxo = -Vxo./Vn;
        Vxo(isnan(Vxo)) = 0;
        Vxo = Vxo.*domDef;
        Vyo = -Vyo./Vn;
        Vyo(isnan(Vyo)) = 0;
        Vyo = Vyo.*domDef;
    otherwise 
        Vxo = zeros(size(X));
        Vyo = zeros(size(X));
end
end

