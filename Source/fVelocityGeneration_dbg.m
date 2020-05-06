function [Vxo, Vyo, PhiC] = fVelocityGeneration_dbg(X, Y, domAt, domBd, domSrc)
%FVELOCITYGENERATION  generates a quiver vector field based on the domains
%of definition, the grid and type
global Dx typeVel

domDef = 1-domBd+domSrc;

switch typeVel
    case 'leftright0'
        Vxo = (X<0)-(X>0);
        Vyo = zeros(size(X));
    case 'centred'
        Vxo = -X./sqrt(X.^2+Y.^2);
        Vxo(isnan(Vxo)) = 0;
        Vyo = -Y./sqrt(X.^2+Y.^2);
        Vyo(isnan(Vyo)) = 0;
    case 'adv-up'
        Vxo = zeros(size(X));
        Vyo = ones(size(X));
    case 'adv-right'
        Vxo = ones(size(X));
        Vyo = zeros(size(X));
    case {'adv-src', 'att-src'}
        Distance = fEikonalCost(X,Y,ones(size(X)),domBd-domAt-domSrc,domAt);
        [Vxo,Vyo] = fDiffFlex(Distance,domDef,Dx,Dx);
        Vn = sqrt(Vxo.^2+Vyo.^2+1-domDef);
        Vxo = -Vxo./Vn;
        Vxo(isnan(Vxo)) = 0;
        Vxo = Vxo.*(domDef);
        Vyo = -Vyo./Vn;
        Vyo(isnan(Vyo)) = 0;
        Vyo = Vyo.*(domDef);
    case {'att','att-adv', 'att2'}
        PhiC = fEikonalCost(X,Y,ones(size(X)),domBd-domAt,domAt);
%         PhiC = X;
        [Vxo,Vyo] = fDiffFlex(PhiC,domDef,Dx,Dx);
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

