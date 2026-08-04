function results = run_nonlinear_test(mpcobj, sys, p, options, L, ...
    idxMV, idxMD, mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, ...
    x0, u0, x_hat0, r_traj, d_phys_traj, Ts)

nu = numel(idxMV);
ny = size(sys.C,1);
Tsim = size(r_traj,1);

x_nl = x0;
u_full = zeros(1, numel(u0));
u_full(:) = u0;   % seed all channels at the validated nominal input
x_hat = x_hat0;

y_log = zeros(Tsim, ny);
u_log = zeros(Tsim, nu);
mpcstate_obj = mpcstate(mpcobj);
mpcstate_obj.Plant = x_hat;

for k = 1:Tsim
    if ~isempty(d_phys_traj)
        u_full(idxMD) = d_phys_traj(k);
    end
    md_norm = ((u_full(idxMD) - mu_u(idxMD)) ./ sigma_u(idxMD)).';

    y_phys = sys_output_wrapper(x_nl, u_full.', p).';
    y_norm = (y_phys - mu_y) ./ sigma_y;
    x_hat = x_hat + L * (y_norm.' - sys.C*x_hat);
    mpcstate_obj.Plant = x_hat;

    if isempty(d_phys_traj)
        mv_norm = mpcmove(mpcobj, mpcstate_obj, y_norm, r_traj(k,:));
    else
        mv_norm = mpcmove(mpcobj, mpcstate_obj, y_norm, r_traj(k,:), md_norm);
    end
    mv_norm = min(max(mv_norm(:), mv_bounds_norm(:,1)), mv_bounds_norm(:,2));
    mv_phys = mv_norm(:).' .* sigma_u(idxMV) + mu_u(idxMV);
    u_full(idxMV) = mv_phys;

    [~, x_traj] = ode15s(@(tt,xx) ode_PEMFC(tt,xx,u_full.'), [0 Ts], x_nl, options);
    if any(isnan(x_traj(:))) || any(isinf(x_traj(:)))
        error('ODE integration failed at k=%d', k);
    end
    x_nl = x_traj(end,:).';
    x_hat = sys.A*x_hat + sys.B(:,idxMV)*mv_norm(:) + sys.B(:,idxMD)*md_norm;

    y_log(k,:) = y_phys;
    u_log(k,:) = mv_phys;
end

results = struct('y', y_log, 'u', u_log, 'x_hat_final', x_hat, 'x_nl_final', x_nl);
end