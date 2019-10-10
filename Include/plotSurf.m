function plotSurf(X,Y,P,Vx,Vy,Axis)
% V1.1 Colormap winter
%% Plot        
    % Faster graph
    contourf(X,Y,min(P,1),40,'EdgeColor','none')
    
    % Higher resolution + 3D
%     surf(X,Y,P,'EdgeColor','none')
%     view(2)

    % Vector field
    % plotQuiver(X,Y,Vx,Vy)
    axis(Axis)
    axis normal
    axis equal
    colorbar
    colormap(winter)
    caxis([0 1])
end
