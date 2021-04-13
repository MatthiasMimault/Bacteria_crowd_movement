function O = fObstacleGeneration(Domain,Attractant,R,r)
%% Generation d'obstacles
% Simple
% O = [-0.1,0.1,-0.1,0.1;0.3,1.5,-0.2,0.1]

% Shifted Cartesian quadrant: 
% xo = Domain(1)+1.3*r:R:Domain(2)-1.1*r;
% yo = Domain(3)+1.3*r:R:Domain(4)-1.1*r;

% True hexagonal lattice
H = sqrt(3)/2*R;
xo = Domain(1)+1.3*r:R:Domain(2)-1.3*r;
yo = Domain(3)+1.3*r:H:Domain(4)-1.3*r;

O = zeros(length(xo)*length(yo),4);
shift = 0;
for j = 1:length(yo)
    for i = 1:length(xo)
        O((j-1)*length(xo)+i,:)=...
            [xo(i)-r+shift*R/2,xo(i)+r+shift*R/2,yo(j)-r,yo(j)+r];
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
    if (C(1)<Domain(1)+1.1*r || C(1)>Domain(2)-1.1*r ...
            || C(2)<Domain(3)+1.1*r || C(2)>Domain(4)-1.1*r)||...
            (C(1)>Attractant(1)-1.1*r && C(1)<Attractant(2)+1.1*r ...
            && C(2)>Attractant(3)-1.25*r && C(2)<Attractant(4)+1.1*r)
        O(n,:)=[];
    end
end
end