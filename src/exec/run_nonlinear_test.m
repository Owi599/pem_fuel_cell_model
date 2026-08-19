function results = run_nonlinear_test(mpcobj, sys, p, options, ...
    idxMV, idxMD, mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, ...
    x0, u0, r_traj, d_phys_traj, Ts, idxUD_dist, ud_phys_traj)

if nargin < 17 || isempty(idxUD_dist)
    idxUD_dist = [];
end
if nargin < 18 || isempty(ud_phys_traj)
    ud_phys_traj = [];
end

nu = numel(idxMV);
ny = size(sys.C, 1);
nx_nl = numel(x0);          % <-- add: full nonlinear state dimension
Tsim = size(r_traj, 1);

x_nl = x0;
u_full = u0(:).';

mpcstate_obj = mpcstate(mpcobj);

y_log = zeros(Tsim, ny);
u_log = zeros(Tsim, nu);
x_log = zeros(Tsim, nx_nl);   % <-- add: full nonlinear state trajectory

nx_mpc = numel(mpcstate_obj.Plant);
xhat_log = nan(Tsim, nx_mpc);

for k = 1:Tsim
    %% 1. Set measured disturbance for both nonlinear plant and MPC
    if isempty(d_phys_traj)
        md_phys = u0(idxMD);
    else
        md_phys = d_phys_traj(k,:);
    end
    u_full(idxMD) = md_phys;

    %% 2. Apply S5 plant-only unmeasured disturbance
    if ~isempty(idxUD_dist)
        u_full(idxUD_dist) = ud_phys_traj(k,:);
    end

    %% 3. Measure nonlinear PEMFC output before applying the new MV
    y_phys = sys_output_wrapper(x_nl, u_full.', p).';
    y_norm = (y_phys - mu_y) ./ sigma_y;

    %% 4. Normalize measured disturbance for MPC
    md_norm = (md_phys - mu_u(idxMD)) ./ sigma_u(idxMD);
    md_norm = reshape(md_norm, 1, []);

    %% 5. MPC move calculation with its internal estimator
    if isempty(idxMD)
        mv_norm = mpcmove(mpcobj, mpcstate_obj, y_norm, r_traj(k,:));
    else
        mv_norm = mpcmove(mpcobj, mpcstate_obj, y_norm, r_traj(k,:), md_norm);
    end

    %% 6. Final numerical saturation guard
    mv_norm = min(max(mv_norm(:), mv_bounds_norm(:,1)), mv_bounds_norm(:,2));

    %% 7. Convert normalized MV command to physical values
    mv_phys = mv_norm.' .* sigma_u(idxMV) + mu_u(idxMV);
    u_full(idxMV) = mv_phys;

    %% 8. Propagate nonlinear DAE plant one sample
    [~, x_traj] = ode15s(@(tt, xx) ode_PEMFC(tt, xx, u_full.'), [0 Ts], x_nl, options);
    if any(~isfinite(x_traj(:)))
        error('ODE integration failed at sample k = %d.', k);
    end
    x_nl = x_traj(end,:).';

    %% 9. Store physical outputs, physical MVs, MPC state estimate, and full state
    y_log(k,:) = y_phys;
    u_log(k,:) = mv_phys;
    xhat_log(k,:) = mpcstate_obj.Plant(:).';
    x_log(k,:) = x_nl.';        % <-- add
end

results = struct( ...
    'y', y_log, ...
    'u', u_log, ...
    'x', x_log, ...             % <-- add: full nonlinear state trajectory, Tsim x nx_nl
    'x_hat_final', mpcstate_obj.Plant, ...
    'x_hat_log', xhat_log, ...
    'x_nl_final', x_nl);
end