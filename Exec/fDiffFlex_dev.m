% Autonatic differentiation of a cartesian vector field, given a definition
% domain. Choices between central, backward and forward diff.
%
%MIMAULT Matthias
%2020
function [Vx,Vy] = fDiffFlex_dev(U,D,dx,dy)
% Definition Boundary
B = 1-D;
[nx,ny] = size(B);

BL = [B(:,1:end-1).*D(:,2:end),zeros(ny,1)];
BR = [zeros(ny,1),B(:,2:end).*D(:,1:end-1)];
UL = [U(:,2:end),zeros(ny,1)].*BL;
UR = [zeros(ny,1),U(:,1:end-1)].*BR;
Upx = [U.*D+UR,zeros(ny,1)];
Umx = [zeros(ny,1),U.*D+UL];

Vc_x = (Upx(:,2:end) - U)/dx;

B;