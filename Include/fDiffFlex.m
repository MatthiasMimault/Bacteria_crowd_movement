%Differentiation of a vectorfield, given an obstace+walls matrix Uobs
%Different choices between central, backward and forward diff.
%
%MIMAULT Matthias
%2014

function [Vx,Vy] = fDiffFlex(V,C,dx,dy)
%V is the value to be differentiate
%C the region of definition of this value
% C=1.*(Gamma>0);
% 
% Bx=([C(:,2:end),zeros(W,1)]-C)>0;
% Fx=([zeros(W,1),C(:,1:end-1)]-C)>0;
% By=([C(2:end,:);zeros(1,L)]-C)>0;
% Fy=([zeros(1,L);C(1:end-1,:)]-C)>0;

Bx=(C-[C(:,2:end),C(:,end)])>0;
Fx=(C-[C(:,1),C(:,1:end-1)])>0;
By=(C-[C(2:end,:);C(end,:)])>0;
Fy=(C-[C(1,:);C(1:end-1,:)])>0;


%Central differentiation
Vc_x = ([V(:,2:end) V(:,end)]-[V(:,1) V(:,1:end-1)])/2/dx;
Vc_y = ([V(2:end,:);V(end,:)]-[V(1,:);V(1:end-1,:)])/2/dy;
% figure
% quiver(Vc_x./(abs(Vc_x)+abs(Vc_y)),Vc_y./(abs(Vc_x)+abs(Vc_y)))

%Forward differentiation
Vf_x = ([V(:,2:end) V(:,end)]-V)/dx;
Vf_y = ([V(2:end,:);V(end,:)]-V)/dy;
% figure
% quiver(Vf_x./(abs(Vf_x)+abs(Vf_y)),Vf_y./(abs(Vf_x)+abs(Vf_y)))

%Backward differentiation
Vb_x = (V-[V(:,1) V(:,1:end-1)])/dx;
Vb_y = (V-[V(1,:);V(1:end-1,:)])/dy;
% figure
% quiver(Vb_x./(abs(Vb_x)+abs(Vb_y)),Vb_y./(abs(Vb_x)+abs(Vb_y)))
% 
% figure
% contourf(C)
% figure
% contourf(C-Bx-Fx)
% figure
% contourf(Fx)
% figure
% contourf(Bx)

Vx = Vc_x.*(C-Bx-Fx)+Vf_x.*Fx+Vb_x.*Bx;
Vy = Vc_y.*(C-By-Fy)+Vf_y.*Fy+Vb_y.*By;
% figure
% quiver(Vx./(abs(Vx)+abs(Vy)),Vy./(abs(Vx)+abs(Vy)))
