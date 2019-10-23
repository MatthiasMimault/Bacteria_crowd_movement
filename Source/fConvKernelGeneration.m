function [E,Ex,Ey] = fConvKernelGeneration()
% fConvolutionkernels produces a normalised centered kernel with finite
% definition based on a radius and a cartesian discretisatrion
global R Dx Dy

[Xe,Ye] = meshgrid(-R:Dx:R, -R:Dy:R);
eta = @(x,y) (1-(x/R).^2).^3.*(1-(y/R).^2).^3;
etax = @(x,y) 3*(1-(x/R).^2).^2.*(1-(y/R).^2).^3.*(-2*x/R.^2);
etay = @(x,y) 3*(1-(x/R).^2).^3.*(1-(y/R).^2).^2.*(-2*y/R.^2);

E = eta(Xe,Ye).*(Xe>-R).*(Xe<R).*(Ye>-R).*(Ye<R);
SE = sum(sum(E));
E=E./SE;
Ex = etax(Xe,Ye).*(Xe>-R).*(Xe<R).*(Ye>-R).*(Ye<R);
Ex=Ex./SE;
Ey = etay(Xe,Ye).*(Xe>-R).*(Xe<R).*(Ye>-R).*(Ye<R);
Ey=Ey./SE;
end

