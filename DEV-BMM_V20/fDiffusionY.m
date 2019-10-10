function b = fDiffusionY(b,Dt,Dx,C,PhiDef,PhiBd)
% Compute the diffusion of b with coefficient C
%% Parameters
s = size(b);
nx = s(1);

% bm = [zeros(s(1),1) b(:,1:end-1)];
% bp = [b(:,2:end) zeros(s(1),1)];
% bm = [b(:,1) b(:,1:end-1)];
% bp = [b(:,2:end) b(:,end)];

%% Neumann boundary
% contourf((PhiBd+[zeros(ny,1) PhiDef(:,1:end-1)]).*[zeros(ny,1) b(:,1:end-1)])
% bp = [zeros(ny,1) PhiDef(:,1:end-1)].*PhiBd.*[zeros(ny,1) b(:,1:end-1)];
% bm = [PhiDef(:,2:end) zeros(ny,1)].*PhiBd.*[b(:,2:end) zeros(ny,1)];
bu = [zeros(1,nx);b(1:end-1,:)].*PhiBd+b;
bd = [b(2:end,:);zeros(1,nx)].*PhiBd+b;


b = b+Dt/Dx/Dx*C*(bd-2*b+bu);
b = b.*PhiDef;
end

