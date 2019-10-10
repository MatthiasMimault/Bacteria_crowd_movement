function b = fInitDatGeneration2(InitSpace,InitValue,X,Y,PhiDef,type)
%UNTITLED Generation of a square initial data defined on PhiDef

% Grid parameters
Dx = X(1,2)-X(1,1);
Dy = Y(2,1)-Y(1,1);
sigma = 0.5;
mu = InitValue;

% Initial data generation
switch type
    case 'stochastic'
        s = size(X);
        v = sigma*randn(s)+mu;     
        b = v.*PhiDef.*(X>InitSpace(1)).*(X<InitSpace(2))...
            .*(Y>InitSpace(3)).*(Y<InitSpace(4));
        b = b.*(b>0).*(b<1);
    otherwise 
        b = InitValue*PhiDef.*(X>InitSpace(1)).*(X<InitSpace(2))...
            .*(Y>InitSpace(3)).*(Y<InitSpace(4));
end
end

