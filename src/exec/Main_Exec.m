%% Main Exec file
close all; clear; clc;
ensure_ode_pemfc_fresh();
%% Load dataset
inputFileName = 'dataset_100_251117_v06';
initial_input = [2.1, 560, 30, 70, 2.3, 164.73, 30, 70, 70, 0, 426, 2.7, 2.4];
Ts = 1;

%% Linearize
[sys, x0, x0_est, u0, mu_u, sigma_u, mu_y, sigma_y, U_Normalized, ...
 y_Normalized, t, idxMV, idxMD, idxUD, p, U] = ...
 linearize_pemfc_N4SID(inputFileName, initial_input);

A = sys.A; B = sys.B; C = sys.C; D = sys.D;
nx_id = size(A,1);

outNames = {'$T_S(z = 123.1 mm) [K]$','$\overline{a}_{H2O} [-]$'};   
varNames = p.inputs.variableNames;

%% Open-loop comparison
y_lin = lsim(sys, U_Normalized, t, x0_est);
figure('Name','Linear (N4SID) vs Nonlinear (DAE) open-loop response');
for kk = 1:2
    
    subplot(2,1,kk);
    plot(t, y_Normalized(:,kk), 'g', 'LineWidth', 2); hold on;grid on;
    plot(t, y_lin(:,kk), 'r--', 'LineWidth', 2);
    ax = gca;
    ax.FontSize= 20;
    ylabel(outNames{kk},'FontSize',24,'FontWeight','bold','Interpreter','latex'); 
    xlabel('Time [s]','FontSize',24,'FontWeight','bold','Interpreter','latex');
    legend('Nonlinear DAE','Linear (N4SID)','FontSize',18,'FontWeight','bold','Location','best');
end
data = iddata(y_Normalized, U_Normalized, Ts);
figure('Name','compare(): NRMSE fit');
compare(data, sys);
%% Role assignment printout
fprintf('\nMV: %s\n', strjoin(varNames(idxMV), ', '));
fprintf('MD: %s\n', strjoin(varNames(idxMD), ', '));
fprintf('Fixed: %s\n', strjoin(varNames(idxUD), ', '));
nu = numel(idxMV); nd = numel(idxMD); ny = 2;
Bu = B(:, idxMV); Bd = B(:, idxMD);

%% Design MPC (default estimator — needed for linear sim tests)
mv_bounds = [ 81.0  534.7;
               1.03  127.0;
              26.9  235.9;
               0.92  139.5;
             305.0  352.7];
[mpcobj, mv_bounds_norm, range_u_norm] = design_mpc(sys, idxMV, idxMD, idxUD, mv_bounds, mu_u, sigma_u, Ts);

%% Setpoint tracking on linear system
Tsim1 = 200;
r1 = y_Normalized(1,:) + [0.02*sigma_y(1), 2*sigma_y(2)];
r1_traj = repmat(r1, Tsim1, 1);
simopt = mpcsimopt(mpcobj);
simopt.PlantInitialState = x0_est;
[y1, t1, u1] = sim(mpcobj, Tsim1, r1_traj, [], simopt);
figure('Name','Test 1: Setpoint Tracking');
for kk = 1:2
    subplot(2,1,kk);
    plot(t1, y1(:,kk), 'b','LineWidth',1.2); hold on; grid on;
    yline(r1(kk),'k--');
    ylabel(outNames{kk}, 'Interpreter','latex'); xlabel('Time [s]');
    legend('Output','Reference');
end

%% Disturbance rejection on linear system
Tsim2 = 200;
r2 = y_Normalized(1,:);
r2_traj = repmat(r2, Tsim2, 1);
d2_traj = zeros(Tsim2,1);
d2_traj(5:20) = 2; d2_traj(20:100) = 3; d2_traj(100:200) = 1;
[y2s, t2, u2s] = sim(mpcobj, Tsim2, r2_traj, d2_traj);
figure('Name','Test 2: Disturbance Rejection');
for kk = 1:2
    subplot(2,1,kk);
    plot(t2, y2s(:,kk), 'r','LineWidth',1.2); hold on; grid on;
    yline(r2(kk),'g--');
    ylabel(outNames{kk}, 'Interpreter','latex'); xlabel('Time [s]');
    legend('Output','Nominal Reference');
end

%% Switch to custom estimator for nonlinear tests and build Kalman gains
setEstimator(mpcobj, 'custom');
Qkf = 1e-4*eye(nx_id);
Rkf = 1e-4*eye(ny);
L = build_kalman_gain(sys, idxMV, idxMD, Qkf, Rkf);

options.Mass = p.M();
options.RelTol = 1e-6;
options.AbsTol = 1e-6;
options.MStateDependence = 'none';

%% Warm-up phase (open loop, settles x_hat)
mpcstate_obj = mpcstate(mpcobj);
x_hat = x0_est;
x_nl = x0;
u_full = zeros(1, size(U,2));
u_full(:) = u0;

warmup_steps = 30;
for k = 1:warmup_steps
    y_phys = sys_output_wrapper(x_nl, u_full.', p).';
    y_norm = (y_phys - mu_y) ./ sigma_y;
    x_hat  = x_hat + L * (y_norm.' - C*x_hat);
    [~, x_traj] = ode15s(@(tt,xx) ode_PEMFC(tt,xx,u_full.'), [0 Ts], x_nl, options);
    x_nl = x_traj(end,:).';
    mv0_norm = ((u_full(idxMV) - mu_u(idxMV)) ./ sigma_u(idxMV)).';
    md0_norm = ((u_full(idxMD) - mu_u(idxMD)) ./ sigma_u(idxMD)).';
    x_hat = A*x_hat + Bu*mv0_norm + Bd*md0_norm;
end
x_hat_warm = x_hat;   % <-- save this, needed by run_nonlinear_test below

%% Closed-loop nominal test (base case)
Tsim_nl = 500;
r_nl_traj = repmat(y_Normalized(1,:), Tsim_nl, 1);
resBase = run_nonlinear_test(mpcobj, sys, p, options, L, idxMV, idxMD, ...
    mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, x0, u0, x_hat_warm, ...
    r_nl_traj, [], Ts);

%% Scenario A: Setpoint tracking on nonlinear plant
r_sp_traj = repmat(y_Normalized(1,:), 300, 1);
r_sp_traj(50:end,1) = r_sp_traj(1,1) + 0.5;
resA = run_nonlinear_test(mpcobj, sys, p, options, L, idxMV, idxMD, ...
    mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, x0, u0, x_hat_warm, ...
    r_sp_traj, [], Ts);

r_sp_phys = r_sp_traj .* sigma_y + mu_y;
t_sp = (0:299)'*Ts;
plot_results(t_sp, resA.y, r_sp_phys, resA.u, mv_bounds, outNames, varNames(idxMV), ...
    'Scenario A: Setpoint Tracking (Nonlinear Plant)');

%% Scenario B: Disturbance rejection on nonlinear plant (data-scaled MD)
r_d_traj = repmat(y_Normalized(1,:), 300, 1);
md_nom = u0(idxMD);

md_min = prctile(U(:,idxMD), 5);
md_max = prctile(U(:,idxMD), 95);
step_amp = 0.05 * (md_max - md_min);

% Clip so the disturbed value never leaves the [md_min, md_max] envelope
up_val   = min(md_nom + step_amp, md_max);
down_val = max(md_nom - step_amp, md_min);

d_phys_traj = repmat(md_nom, 300, 1);
ramp_len = 10;
up_ramp   = linspace(0, up_val - md_nom, ramp_len)';
down_ramp = linspace(0, md_nom - down_val, ramp_len)';

d_phys_traj(50:59)   = md_nom + up_ramp;
d_phys_traj(60:119)  = up_val;
d_phys_traj(120:129) = up_val - up_ramp;
d_phys_traj(130:199) = down_val;
d_phys_traj(200:209) = down_val + down_ramp;
d_phys_traj(210:end) = md_nom;

resB = run_nonlinear_test(mpcobj, sys, p, options, L, idxMV, idxMD, ...
    mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, x0, u0, x_hat_warm, ...
    r_d_traj, d_phys_traj, Ts);

r_d_phys = r_d_traj .* sigma_y + mu_y;
t_d = (0:299)'*Ts;
plot_results(t_d, resB.y, r_d_phys, resB.u, mv_bounds, outNames, varNames(idxMV), ...
    'Scenario B: Disturbance Rejection (Nonlinear Plant)');