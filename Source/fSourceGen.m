function PhiSrc = fSourceGen(X,Y,PhiAt,Space,Attractant,typeBdy)
%% Initialisation
Dx = X(1,2)-X(1,1);
Dy = Y(2,1)-Y(1,1);

%% Source bdy design
switch typeBdy
    case 'src'
        PhiSrc = ones(size(X))...
        -(X>Space(1)-Dx/2).*(X<Space(2)+Dx/2)...
        .*(Y>Space(3)-Dy/2).*(Y<Space(4)+Dy/2);
        PhiSrc = PhiSrc.*(Y<0);
    
    otherwise
        PhiSrc = zeros(size(X));
end
        