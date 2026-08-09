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

b_u = U(1,:);
b_y = y(1,:);

sc_u= [1000, 1000, 1000, 1000, 1000, 100, 100, 100, ...
        100000, 100000, 1];
sc_y= [100, 1];

assert(numel(sc_u) == size(U,2), ...
    'sc_u must contain exactly one scale factor for every input.');
assert(numel(sc_y) == size(y,2), ...
    'sc_y must contain exactly one scale factor for every output.');
assert(all(sc_u > 0) && all(sc_y > 0), ...
    'All normalization scale factors must be strictly positive.');

U_n = (U - b_u) ./ sc_u;
y_n = (y - b_y) ./ sc_y;

mu_u = b_u;
sigma_u = sc_u;
mu_y = b_y;
sigma_y = sc_y;

Ts = mean(diff(t));
assert(all(abs(diff(t) - Ts) < 1e-10), ...
    'N4SID requires a uniformly sampled dataset.');

opt = n4sidOptions('N4Weight','auto','Focus','simulation', ...
    'N4Horizon',[20 40 40], 'EnforceStability',true, 'InitialState','zero','OutputWeight',diag([0.3,1]));
[sys, x0_est] = n4sid(U_n, y_n, 6, 'Ts', Ts, opt);
sys.UserData.x0_est = x0_est;   % stash for later reuse instead of a separate output
sys.UserData.t = t; sys.UserData.U_n = U_n; sys.UserData.y_n = y_n;  % for diagnostics/plots

idxMV = [1 2 3 4 8];
idxMD = 11;
idxUD = setdiff(1:size(U,2), [idxMV idxMD]);

fprintf('Stable: %d | Controllable: %d | Observable: %d | Order: %d\n', ...
    all(abs(eig(sys.A))<1), rank(ctrb(sys.A,sys.B))==size(sys.A,1), ...
    rank(obsv(sys.A,sys.C))==size(sys.A,1), size(sys.A,1));

end