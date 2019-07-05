function plotBd(X,Y,P,Axis)
% V1.1 Colormap winter
%% Plot        
    % Faster graph
    contour(X,Y,P,'EdgeColor','magenta')
    
    % Higher resolution + 3D
%     surf(X,Y,P,'EdgeColor','none')
%     view(2)

    % Vector field
    % plotQuiver(X,Y,Vx,Vy)
%     axis(Axis)
%     axis normal
%     axis equal
end

