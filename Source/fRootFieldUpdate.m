function [Vx,Vy] = fRootFieldUpdate(Vxo,Vyo,typeRepulsion)
if typeRepulsion == 1
    Vx = -Vxo;
    Vy = -Vyo;
else
    Vx = Vxo;
    Vy = Vyo;
end

