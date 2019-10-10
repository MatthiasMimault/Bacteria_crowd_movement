function b = fDiffusionX(b,Dt,Dx,C,PhiDef,PhiBd)
% Compute the diffusion of b with coefficient C

%% Parameters
s = size(b);
ny = s(2);

% bm = [zeros(s(1),1) b(:,1:end-1)];
% bp = [b(:,2:end) zeros(s(1),1)];
% bm = [b(:,1) b(:,1:end-1)];
% bp = [b(:,2:end) b(:,end)];

%% Neumann boundary
% contourf((PhiBd+[zeros(ny,1) PhiDef(:,1:end-1)]).*[zeros(ny,1) b(:,1:end-1)])
% bp = [zeros(ny,1) PhiDef(:,1:end-1)].*PhiBd.*[zeros(ny,1) b(:,1:end-1)];
% bm = [PhiDef(:,2:end) zeros(ny,1)].*PhiBd.*[b(:,2:end) zeros(ny,1)];
bp = [zeros(ny,1) b(:,1:end-1)].*PhiBd+b;
% contourf(bp)
% pause
bm = [b(:,2:end) zeros(ny,1)].*PhiBd+b;




b = b+Dt/Dx/Dx*C*(bm-2*b+bp);
b = b.*PhiDef;
end

