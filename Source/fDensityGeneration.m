function b = fDensityGeneration(X,Y,value,dom)
%FDENSITYGENERATION generates the density b according to a cartesian region
%domInit, an initial value initValue and a type typeInit

global typeInit A C

% Stochastic parameters
sigma = 0.5;
mu = value;

switch typeInit
    case 'stochastic'
        s = size(X);
        v = sigma*randn(s)+mu;     
        b = v.*(X>dom(1)).*(X<dom(2)).*(Y>dom(3)).*(Y<dom(4));
        b = b.*(b>0).*(b<1);
    case 'square'
        b = value*(X>dom(1)).*(X<dom(2)).*(Y>dom(3)).*(Y<dom(4));
    case 'solstat'
        b = exp(-A/C*sqrt(X.^2+Y.^2));
    case 'solstat1D'
        b = exp(-A/C*abs(X));
    case 'solstat2D'
        b = exp(-A/C*sqrt(X.^2+Y.^2));
    otherwise 
        b = zeros(size(X));
end
end

