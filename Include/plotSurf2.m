function plotSurf2(X,Y,P,Vx,Vy,Axis)
% V1.1 Colormap winter
%% Plot        
    Pmax = 20000;
    % Faster graph
    P(P>Pmax) = Pmax;
    P(P<0.1) = NaN;
%     contourf(X,Y,P,logspace(-1,5,100),'EdgeColor','none')
    contourf(X,Y,P,logspace(log10(0.1),log10(Pmax),100),'EdgeColor','none')
    
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
    caxis([0.1 Pmax])
    set(gca,'colorscale','log')
end

