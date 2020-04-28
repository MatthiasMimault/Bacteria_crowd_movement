function plotSurf2(X,Y,P,Vx,Vy,Axis)
% V1.1 Colormap winter
%% Plot        
    Pmin = 0;
    Pmax = 200;
    % Faster graph
    P(P>Pmax) = Pmax;
    P(P<=Pmin) = NaN;
    contourf(X,Y,P,100,'EdgeColor','none')
%     contourf(X,Y,P,logspace(log10(Pmin),log10(Pmax),100),'EdgeColor','none')
    
    % Higher resolution + 3D
%     surf(X,Y,P,'EdgeColor','none')
%     view(2)

    % Vector field
    % plotQuiver(X,Y,Vx,Vy)
    axis(Axis)
    axis normal
    axis equal
    colorbar
    colormap(jet)
    caxis([Pmin Pmax])
%     set(gca,'colorscale','log')
end

