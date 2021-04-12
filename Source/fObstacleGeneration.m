function O = fObstacleGeneration(Domain,Attractant,R,r)
%% Generation d'obstacles
% Simple
% O = [-0.1,0.1,-0.1,0.1;0.3,1.5,-0.2,0.1]

% Serial Cartesian quadrant: 
% R = 0.5; r = 0.2;
%%%xo = Domain(1)+R-r/2:R:Domain(2)-R+r/2;
%%%yo = Domain(3)+R-r/2:R:Domain(4)-R+r/2;
xo = Domain(1):R:Domain(2);
yo = Domain(3):R:Domain(4);
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

%% Cleaning
% dev = 3;
s = size(O);
% RadiusObstacle = 
for n = s(1):-1:1
    % Center and Radius of obstacle n
    C = [(O(n,1)+O(n,2))/2, (O(n,3)+O(n,4))/2];
    if (C(1)<Domain(1)+2*r || C(1)>Domain(2)-2*r ...
            || C(2)<Domain(3)+2*r || C(2)>Domain(4)-2*r)||...
            (C(1)>Attractant(1)-2*r && C(1)<Attractant(2)+2*r ...
            && C(2)>Attractant(3)-2*r && C(2)<Attractant(4)+2*r)
        O(n,:)=[];
    end
end
end