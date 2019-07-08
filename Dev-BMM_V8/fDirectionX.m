function Vx = fDirectionX(b,Vxo,Vyo,Ex,Ey,PhiDef)

%% Parameters
epsilon = 0.5;

%% Vx computation
% Perceived density
UEx = convolve2(b.*PhiDef,Ex,'same');
UEy = convolve2(b.*PhiDef,Ey,'same');
SU = sqrt(1+UEx.^2+UEy.^2);

% Deviation
Ix = epsilon*UEx./SU;
Iy = epsilon*UEy./SU;

% Vx update
Vn = sqrt((Vxo+Ix).^2+(Vyo+Iy).^2);
Vx = (Vxo+Ix)./Vn;
Vy = (Vyo+Iy)./Vn;
Vx(isnan(Vx)) = 0;
Vx(isnan(Vy)) = 0;
Vx = Vx.*PhiDef;
end

