function plotQuiver(X,Y,Vx,Vy)

%% 0a. Prelude
% hold on
% if exist(['Data/' namePathTest '/DensityPNG/'], 'dir')==0
%     mkdir(['Data/' namePathTest '/DensityPNG/']);
% end
% if exist(['Data/' namePathTest '/DensityFIGURES/'], 'dir')==0
%     mkdir(['Data/' namePathTest '/DensityFIGURES/']);
% end

%% 0b. Ratio d'affichage
[Lx,Ly]=size(X);
%Plus c est haut, plus c est ressere
K=320;
kx = ceil(Lx/K);
ky = ceil(Ly/K);

%% 1. Plot Quiver (Grey)
quiver(X(1:kx:end,1:ky:end),Y(1:kx:end,1:ky:end),...
    Vx(1:kx:end,1:ky:end),Vy(1:kx:end,1:ky:end),'color',[0.5 0.5 0.5])
axis equal
%% 3. Save figure
% saveas(gcf,['Data/' namePathTest '/DensityPNG/Quiver' int2str(n)],'png')
% savefig(['Data/' namePathTest '/DensityFIGURES/Quiver' int2str(n)])

end

