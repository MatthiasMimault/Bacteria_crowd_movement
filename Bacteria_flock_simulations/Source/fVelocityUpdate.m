function [Vx, Vy] = fVelocityUpdate(b,Vxo,Vyo,Ex,Ey,Bdy)
%fVelocityUpdate updates the vecto field V with the density estimation
%around the bacteria (attraction)

global typeVel

switch typeVel
    case 'att2'
        [Vx, Vy] = fDeviation_dev(b,Vxo,Vyo,Ex,Ey,Bdy);
    case {'att', 'att-src'}
        [Vx, Vy] = fDeviation(b,Vxo,Vyo,Ex,Ey,Bdy);
    
    otherwise
        Vx = Vxo;
        Vy = Vyo;
end
end

function [Vx, Vy] = fDeviation(b,Vxo,Vyo,Ex,Ey,Bdy)
%% Parameters
global kappa
Def = 1-Bdy;

%% Vx computation
% Perceived density
UEx = convolve2(b.*Def,Ex,'same');
UEy = convolve2(b.*Def,Ey,'same');
SU = sqrt(1+UEx.^2+UEy.^2);

% Deviation
Ix = kappa*UEx./SU.*Def;
Iy = kappa*UEy./SU.*Def;

% V update
Vn = sqrt((Vxo+Ix).^2+(Vyo+Iy).^2);

Vx = (Vxo+Ix)./Vn;
Vy = (Vyo+Iy)./Vn;

% V cleaning of nan values
nanVx = isnan(Vx);
nanVy = isnan(Vy);

Vx(nanVx) = 0;
Vx(nanVy) = 0;
Vx = Vx.*Def;

Vy(nanVx) = 0;
Vy(nanVy) = 0;
Vy = Vy.*Def;
end

function [Vx, Vy] = fDeviation_dev(b,Vxo,Vyo,Ex,Ey,Bdy)
%% Parameters
global kappa
Def = 1-Bdy;

%% Vx computation
% Perceived density
UEx = convolve2(b.*Def,Ex,'same');
UEy = convolve2(b.*Def,Ey,'same');
SU = sqrt(1+UEx.^2+UEy.^2);

% Deviation
Ix = kappa*UEx./SU.*Def;
Iy = kappa*UEy./SU.*Def;

% V update
Vx = Vxo+Ix;
Vy = Vyo+Iy;

Vx = Vx.*Def;
Vy = Vy.*Def;
end