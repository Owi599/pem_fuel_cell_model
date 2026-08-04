function [sys, x0, x0_est, u0, mu_u, sigma_u, mu_y, sigma_y, U_n, y_n , t, idxMV, idxMD, idxUD,p,U] = linearize_pemfc_N4SID(datasetName, in1)

p = mod_param_PEMFC();
options.Mass = p.M();
options.RelTol = 1e-6;
options.AbsTol = 1e-6;
options.MStateDependence = 'none';

[u_traj, ~] = loadMatFile([datasetName,'_inputs.mat'], [], p.inputs.variableNames);
uInterp = griddedInterpolant(u_traj.time, u_traj.data, 'pchip', 'nearest');

initialInput = p.testbench2struct(in1.');
x0_guess = steady_state_PEMFC(p, initialInput, options);

tspan = u_traj.time;
U = u_traj.data;
if size(U,1) ~= numel(tspan); U = U.'; end

u = @(t) uInterp(t).';
[t, x] = ode15s(@(t,x) ode_PEMFC(t,x,u(t)), tspan, x0_guess, options);

x0 = x(1,:).';        % ode15s-validated consistent state
u0 = U(1,:);           % its matching input

N = size(x,1);
y = zeros(N,2);
for k = 1:N
    y(k,:) = sys_output_wrapper(x(k,:).', U(k,:).', p).';
end

mu_u = mean(U,1); sigma_u = std(U,1); sigma_u(sigma_u==0) = 1;
mu_y = mean(y,1); sigma_y = std(y,1); sigma_y(sigma_y==0) = 1;

U_n = (U - mu_u) ./ sigma_u;
y_n = (y - mu_y) ./ sigma_y;

Ts = 1;
opt = n4sidOptions('N4Weight','auto','Focus','prediction', ...
    'N4Horizon',[11 100 100], 'EnforceStability',true, 'InitialState','zero');
[sys, x0_est] = n4sid(U_n, y_n, 10:20, 'Ts', Ts, opt);
sys.UserData.x0_est = x0_est;   % stash for later reuse instead of a separate output
sys.UserData.t = t; sys.UserData.U_n = U_n; sys.UserData.y_n = y_n;  % for diagnostics/plots

idxMV = [1 2 3 4 8];
idxMD = 11;
idxUD = setdiff(1:size(U,2), [idxMV idxMD]);

fprintf('Stable: %d | Controllable: %d | Observable: %d | Order: %d\n', ...
    all(abs(eig(sys.A))<1), rank(ctrb(sys.A,sys.B))==size(sys.A,1), ...
    rank(obsv(sys.A,sys.C))==size(sys.A,1), size(sys.A,1));

end