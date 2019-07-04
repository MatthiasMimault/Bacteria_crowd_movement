function O = fObstacleGeneration(Domain,R,r)
%% Generation d'obstacles
% Simple
% O = [-0.1,0.1,-0.1,0.1;0.3,1.5,-0.2,0.1]

% Serial Cartesian quadrant: 
% R = 0.5; r = 0.2;
xo = Domain(1)+R/2:R:Domain(2)-R/2;
yo = Domain(3)+R/2:R:Domain(4)-R/2;
% [Xo,Yo] = meshgrid(xo,yo);
O = zeros(length(xo)*length(yo),4);
shift = 0;
for j = 1:length(yo)
    for i = 1:length(xo)
        O((j-1)*length(xo)+i,:)=...
            [xo(i)-r/2+shift*R/2,xo(i)+r/2+shift*R/2,yo(j)-r/2,yo(j)+r/2];
    end
    if shift < 0.5
        shift = 1;
    else
        shift = 0;
    end
end