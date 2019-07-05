function O = fObstacleGeneration(Domain,Attractant,R,r)
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

%% Cleaning
s = size(O);
for n = s(1):-1:1
    % Center and Radius of obstacle n
    C = [(O(n,1)+O(n,2))/2, (O(n,3)+O(n,4))/2];
    R = O(n,2)-O(n,1);
    if C(1) < Domain(1)+R || C(1) > Domain(2)-R || ...
            C(2) < Domain(3)+R || C(2) > Domain(4)-R || ...
            (C(1) > Attractant(1)-R && C(1) < Attractant(2)+R && ...
            C(2) > Attractant(3)-R && C(2) < Attractant(4)+R)
        O(n,:)=[];
    end
end