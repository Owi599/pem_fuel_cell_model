function L = build_kalman_gain(sys, idxMV, idxMD, Qkf, Rkf)

nx = size(sys.A,1); ny = size(sys.C,1);
Bu = sys.B(:,idxMV); Bd = sys.B(:,idxMD);
Du = sys.D(:,idxMV); Dd = sys.D(:,idxMD);

sys_kf = ss(sys.A, [Bu Bd eye(nx)], sys.C, [Du Dd zeros(ny,nx)], sys.Ts);
[~, L, ~] = kalman(sys_kf, Qkf, Rkf, 0, 1:ny, 1:(numel(idxMV)+numel(idxMD)));

end