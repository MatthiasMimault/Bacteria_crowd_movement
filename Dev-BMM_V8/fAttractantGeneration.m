function PhiAt = fAttractantGeneration(X,Y,typeAt)

%% Initialisation
Dx = X(1,2)-X(1,1);
Dy = Y(2,1)-Y(1,1);

%% Assignation
if strcmp(typeAt,'root')
    Attractant = [-0.2 0.1 0.3 1];
    Xo = (Attractant(1)+Attractant(2))/2; Yo = Attractant(3);
    R = (Attractant(2)-Attractant(1))/2;
    PhiAt = ones(size(X))...
    .*(X>Attractant(1)-Dx/2).*(X<Attractant(2)+Dx/2)...
    .*(Y>Attractant(3)-Dy/2).*(Y<Attractant(4)+Dy/2)...
    +ones(size(X)).*(((X-Xo).^2+(Y-Yo).^2)<=R^2).*(Y<=Attractant(3)-Dy/2);
else
    Attractant = [-1, 1, 1, 1.5];
    PhiAt = ones(size(X))...
        .*(X>Attractant(1)-Dx/2).*(X<Attractant(2)+Dx/2)...
        .*(Y>Attractant(3)-Dy/2).*(Y<Attractant(4)+Dy/2);
end