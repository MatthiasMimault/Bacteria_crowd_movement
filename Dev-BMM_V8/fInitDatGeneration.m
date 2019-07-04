function b = fInitDatGeneration(InitSpace,InitValue,X,Y,PhiDef)
%UNTITLED Generation of a square initial data defined on PhiDef

% Grid parameters
Dx = X(1,2)-X(1,1);
Dy = Y(2,1)-Y(1,1);

% Initial data generation
b = InitValue*PhiDef.*(X>InitSpace(1)).*(X<InitSpace(2))...
    .*(Y>InitSpace(3)).*(Y<InitSpace(4));
end

