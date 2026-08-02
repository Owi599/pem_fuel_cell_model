%% PEMFC: N4SID linear model 
close all; clear; clc;

%% Loading params and ode options
p = mod_param_PEMFC();
options.Mass = p.M();
options.RelTol = 1e-6;
options.AbsTol = 1e-6;
options.MStateDependence = 'none';

%% Load dataset 
inputFileName = 'dataset_100_251117_v06';
[u_traj, ~] = loadMatFile([inputFileName,'_inputs.mat'], [], p.inputs.variableNames);
uInterpolant_pp = griddedInterpolant(u_traj.time, u_traj.data, 'pchip', 'nearest');

in1 = [2.1, 560, 30, 70, 2.3, 164.73, 30, 70, 70, 0, 426, 2.7, 2.4];
initialInput = p.testbench2struct(in1.');
x0 = steady_state_PEMFC(p, initialInput, options);
tspan = u_traj.time;
U = u_traj.data;
if size(U,1) ~= numel(tspan); U = U.'; end

%% Simulate nonlinear DAE
u = @(t) uInterpolant_pp(t).';
[t, x] = ode15s(@(t,x) ode_PEMFC(t,x,u(t)), tspan, x0, options);
x0 = x(1,:).';
u0 = U(1,:);
N = size(x,1);
y = zeros(N,2);
for k = 1:N
    y(k,:) = sys_output_wrapper(x(k,:).', U(k,:).', p).';
end

%% 6.1 compute means and std
mu_u = mean(U,1); sigma_u = std(U,1);
mu_y = mean(y,1); sigma_y = std(y,1);
% guard against any zero-variance columns
sigma_u(sigma_u == 0) = 1;
sigma_y(sigma_y == 0) = 1;

normalize_u = @(U) (U - mu_u) ./ sigma_u;
normalize_y = @(Y) (Y - mu_y) ./ sigma_y;
U_Normalized = normalize_u(U);
y_Normalized = normalize_y(y);
fprintf('\nPer-channel training stats (mean, std):\n');
fprintf('  y1 (T_S(10))    : mean=%.3g, std=%.3g\n', mu_y(1), sigma_y(1));
fprintf('  y2 (a_H2O_avg)  : mean=%.3g, std=%.3g\n', mu_y(2), sigma_y(2));

%% Identify linear model (N4SID, raw physical units) 


Ts = 1;
nx = 10:20;
r_forward  = 11;
sy_outputs = 100;
su_outputs = 100;

horizons = [r_forward sy_outputs su_outputs];

data = iddata(y_Normalized,U_Normalized,Ts);

%OPT = ssestOptions( 'Focus','simulation', ...
%  'EnforceStability', false);
%OPT.OutputWeight =diag([0.1,1]);
opt = n4sidOptions('N4Weight','auto','Focus','prediction', ...
    'N4Horizon',horizons, 'EnforceStability',true, ...
    'InitialState', 'zero');
[sys, x0_est] = n4sid(U_Normalized, y_Normalized, nx, 'Ts', Ts, opt);
A = sys.A; B = sys.B; C = sys.C; D = sys.D;
nx_id = size(A,1);
%sys_refined = ssest(data,sys,OPT);  % refining the system

fprintf('Stable:       %d\n', all(abs(eig(A)) < 1));
fprintf('Controllable: %d\n', rank(ctrb(A,B)) == nx_id);
fprintf('Observable:   %d\n', rank(obsv(A,C)) == nx_id);
fprintf('Model order:  %d\n', nx_id);

%% Linear vs nonlinear open-loop comparison

y_lin = lsim(sys, U_Normalized, t,x0_est);
outNames = {'T_S(10) [K]','\overline{a}_{H2O} [-]'};
figure('Name','Linear (N4SID) vs Nonlinear (DAE) open-loop response');

for kk = 1:2
    subplot(2,1,kk);
    plot(t, y_Normalized(:,kk), 'g', 'LineWidth', 1.2); hold on;
    plot(t, y_lin(:,kk), 'r--', 'LineWidth', 1.2);
    ylabel(outNames{kk}); xlabel('Time [s]');
    legend('Nonlinear DAE','Linear (N4SID)');
end
figure('Name','compare(): NRMSE fit, linear model vs nonlinear DAE data');
compare(data, sys);

%% Role assignment
idx_MV = [1 2 3 4 8];
idx_MD = 11;
idx_UD = [5 6 7 9 10];
varNames = p.inputs.variableNames;
fprintf('\nMV: %s\n', strjoin(varNames(idx_MV), ', '));
fprintf('MD: %s\n', strjoin(varNames(idx_MD), ', '));
fprintf('Fixed: %s\n', strjoin(varNames(idx_UD), ', '));

nu = numel(idx_MV); nd = numel(idx_MD); ny = 2;
Bu = B(:, idx_MV); Bd = B(:, idx_MD);
Du = D(:, idx_MV); Dd = D(:, idx_MD);
u_UD_nom = mean(U(:,idx_UD), 1);
u_UD_nom_norm = (mean(U(:,idx_UD),1) - mu_u(idx_UD)) ./ sigma_u(idx_UD);

mv_bounds = [ 81.0  534.7;
               1.03  127.0;
              26.9  235.9;
               0.92  139.5;
             305.0  352.7];
mv_bounds_norm = (mv_bounds - mu_u(idx_MV)') ./ sigma_u(idx_MV)';
range_u = diff(mv_bounds, 1, 2);
range_u_norm = mv_bounds_norm(:,2) - mv_bounds_norm(:,1);
fprintf('Controllable (MVs only): %d\n', rank(ctrb(A,Bu)) == nx_id);

%% MPC – stripped-down first pass ( just to check if y2 is converging towards value near the op)
plant = sys;

% Explicitly mark which inputs are MVs, everything else as UD (ignored)
allInputs = 1:size(B,2);              % 1:11
idx_MV    = [1 2 3 4 8];              % your MVs
idx_UDloc = setdiff(allInputs, idx_MV);

plant = setmpcsignals(plant, ...
    'MV', idx_MV, ...
    'UD', idx_UDloc);                 % no MD for this minimal test

p_horizon = 40;
m_horizon = 10;
mpcobj = mpc(plant, Ts, p_horizon, m_horizon);

% Very simple weights
mpcobj.Weights.OutputVariables          = [1 10];
mpcobj.Weights.ManipulatedVariables     = 0.01 * ones(1, numel(idx_MV));
mpcobj.Weights.ManipulatedVariablesRate = 0.1  * ones(1, numel(idx_MV));

% No manual ScaleFactor tuning at this stage
% (let MPC auto-scale internally)
for i = 1:numel(idx_MV)
    mpcobj.MV(i).Min = mv_bounds_norm(i,1);
    mpcobj.MV(i).Max = mv_bounds_norm(i,2);
end

review(mpcobj);   % just to make sure it's still OK

%Test 1: hold at nominal, then small step in y1
Tsim1 = 2000;
r_nom   = mu_y;                    % [mu_y1, mu_y2]
r1_traj = repmat(r_nom, Tsim1, 1);

[y1, t1, u1] = sim(mpcobj, Tsim1, r1_traj);  % default initial state

figure('Name','Test 1: Small step around nominal');
for kk = 1:2
    subplot(2,1,kk);
    plot(t1, y1(:,kk), 'b','LineWidth',1.2); hold on;
    yline(r_nom(kk),'k--');
    ylabel(outNames{kk}); xlabel('Time [s]');
    legend('Output','Reference');
end
%% MPC object with tuning 
plant= sys;
plant= setmpcsignals(plant,'MV',idx_MV,'MD',idx_MD,'UD',idx_UD);
plant.InputGroup.Unmeasured = idx_UD;  % or exclude them from disturbance estimation
p_horizon = 100;   % Prediction horizon
m_horizon = 5;    % Control horizon
mpcobj = mpc(plant, Ts, p_horizon, m_horizon);

mpcobj.Weights.OutputVariables = [1,1];
mpcobj.Weights.ManipulatedVariables      = 0.01 * ones(1, nu);    % small usage penalty
mpcobj.Weights.ManipulatedVariablesRate  = 0.1 ./ range_u_norm(:)';    % keep these
mpcobj.Weights.ManipulatedVariablesRate = 100 * mpcobj.Weights.ManipulatedVariablesRate;
for i = 1:nu
    mpcobj.MV(i).ScaleFactor = 0.5*range_u_norm(i);
end
mpcobj.OV(1).ScaleFactor = 1;
mpcobj.OV(2).ScaleFactor = 1;
for i = 1:numel(idx_MV)
    mpcobj.MV(i).Min = mv_bounds_norm(i,1);
    mpcobj.MV(i).Max = mv_bounds_norm(i,2);
end
review(mpcobj);             % sanity check
mpcstate_obj = mpcstate(mpcobj);

%% Set point tracking on linear system
Tsim1 = 200;
r1 = y_Normalized(1,:) + [0.05*sigma_y(1) , 2*sigma_y(2) ];
r1_traj = repmat(r1, Tsim1, 1);
simopt = mpcsimopt(mpcobj);
simopt.PlantInitialState = x0_est;     % nominal linear state matching mu_y
simopt.MVSignal = [];                  % leave default unless you want open-loop MV override

[y1, t1, u1] = sim(mpcobj, Tsim1, r1_traj, [], simopt);
figure('Name','Test 1: Setpoint Tracking');
for kk = 1:2
    subplot(2,1,kk);
    plot(t1, y1(:,kk), 'b','LineWidth',1.2); hold on; grid on ;
    yline(r1(kk),'k--');
    ylabel(outNames{kk},'Interpreter','latex'); xlabel('Time [s]','Interpreter','latex');
    legend('Output','Reference');
end

%% Disturbance rejection on linear system
Tsim2 = 200;
r2 = y_Normalized(1,:);
r2_traj = repmat(r2, Tsim2, 1);        % hold at nominal steady output
d2_traj = zeros(Tsim2,1);
d2_traj(5:20) = 2; % step disturbance in I_cell at t=50
d2_traj(20:100)= 3;
d2_traj(100:200) = 1;
[y2s, t2, u2s] = sim(mpcobj, Tsim2, r2_traj, d2_traj);

figure('Name','Test 2: Disturbance Rejection');
for kk = 1:2
    subplot(2,1,kk);
    plot(t2, y2s(:,kk), 'r','LineWidth',1.2); hold on;grid on;
    yline(r2(kk),'g--');
    ylabel(outNames{kk}); xlabel('Time [s]');
    legend('Output','Nominal Reference');
end
%% Kalman Observer
sys_kf = ss(A, [Bu Bd eye(nx_id)], C, [Du Dd zeros(ny,nx_id)], Ts);
known  = 1:(nu+nd);       
Qkf    = 1e-4*eye(nx_id);
Rkf    = 1e-4*eye(ny);
[~, L, ~] = kalman(sys_kf, Qkf, Rkf, 0, 1:ny, known);

%%
setEstimator(mpcobj,'custom');
for i = 1:nu
    mpcobj.MV(i).RateMin = -0.05 * range_u_norm(i);   % cap per-step change
    mpcobj.MV(i).RateMax =  0.05 * range_u_norm(i);
end

mpcstate_obj = mpcstate(mpcobj);
x_hat = x0_est;
x_nl = x0;
u_full = zeros(1, size(U,2));
u_full(idx_UD) = u0(idx_UD);
u_full(idx_MD) = u0(idx_MD);
u_full(idx_MV) = u0(idx_MV);     % <-- ADD THIS: seed MVs at nominal, not zero

warmup_steps = 30;
for k = 1:warmup_steps
    y_phys = sys_output_wrapper(x_nl, u_full.', p).';
    y_norm = (y_phys - mu_y) ./ sigma_y;
    x_hat  = x_hat + L * (y_norm.' - C*x_hat);
    [~, x_traj] = ode15s(@(tt,xx) ode_PEMFC(tt,xx,u_full.'), [0 Ts], x_nl, options);
    x_nl = x_traj(end,:).';
    mv0_norm = ((u_full(idx_MV) - mu_u(idx_MV)) ./ sigma_u(idx_MV)).';
    md0_norm = ((u_full(idx_MD) - mu_u(idx_MD)) ./ sigma_u(idx_MD)).';
    x_hat = A*x_hat + Bu*mv0_norm + Bd*md0_norm;
end
mpcstate_obj.Plant = x_hat;


%% Step 4: Closed-loop test against nonlinear plant
Tsim_nl = 500;
r_nl = y_Normalized(1,:);
r_nl_traj = repmat(r_nl, Tsim_nl, 1);
y_nl_log = zeros(Tsim_nl, ny);
u_nl_log = zeros(Tsim_nl, nu);

for k = 1:Tsim_nl
    y_phys = sys_output_wrapper(x_nl, u_full.', p).';
    y_norm = (y_phys - mu_y) ./ sigma_y;
    x_hat  = x_hat + L * (y_norm.' - C*x_hat);
    mpcstate_obj.Plant = x_hat;

    [mv_norm, info] = mpcmove(mpcobj, mpcstate_obj, y_norm, r_nl_traj(k,:));
    mv_norm = min(max(mv_norm(:), mv_bounds_norm(:,1)), mv_bounds_norm(:,2));
    mv_phys = mv_norm(:).' .* sigma_u(idx_MV) + mu_u(idx_MV);
    u_full(idx_MV) = mv_phys;

    [~, x_traj] = ode15s(@(tt,xx) ode_PEMFC(tt,xx,u_full.'), [0 Ts], x_nl, options);
    x_nl = x_traj(end,:).';

    x_hat = A*x_hat + Bu*mv_norm(:) + Bd*((u_full(idx_MD)-mu_u(idx_MD))./sigma_u(idx_MD)).';
    y_nl_log(k,:) = y_phys;
    u_nl_log(k,:) = mv_phys;
end
%% Scenario A: Setpoint tracking on nonlinear plant
x_hat_A = x_hat;      % branch off from the warmed-up state
x_nl_A  = x_nl;
u_full_A = u_full;

Tsim_sp = 300;
r_sp = y_Normalized(1,:);           % start at nominal
r_sp_traj = repmat(r_sp, Tsim_sp, 1);
r_sp_traj(50:end, 1) = r_sp(1) + 0.5;   % step in y1 (normalized units) at k=50

y_sp_log = zeros(Tsim_sp, ny);
u_sp_log = zeros(Tsim_sp, nu);

for k = 1:Tsim_sp
    y_phys = sys_output_wrapper(x_nl_A, u_full_A.', p).';
    y_norm = (y_phys - mu_y) ./ sigma_y;
    x_hat_A = x_hat_A + L * (y_norm.' - C*x_hat_A);
    mpcstate_obj.Plant = x_hat_A;
    [mv_norm, ~] = mpcmove(mpcobj, mpcstate_obj, y_norm, r_sp_traj(k,:));
    mv_norm = min(max(mv_norm(:), mv_bounds_norm(:,1)), mv_bounds_norm(:,2));
    mv_phys = mv_norm(:).' .* sigma_u(idx_MV) + mu_u(idx_MV);
    u_full_A(idx_MV) = mv_phys;
    [~, x_traj] = ode15s(@(tt,xx) ode_PEMFC(tt,xx,u_full_A.'), [0 Ts], x_nl_A, options);
    x_nl_A = x_traj(end,:).';
    x_hat_A = A*x_hat_A + Bu*mv_norm(:) + Bd*((u_full_A(idx_MD)-mu_u(idx_MD))./sigma_u(idx_MD)).';
    y_sp_log(k,:) = y_phys;
    u_sp_log(k,:) = mv_phys;
end

r_sp_phys = r_sp_traj .* sigma_y + mu_y;
t_sp = (0:Tsim_sp-1)'*Ts;

figure('Name','Scenario A: Setpoint Tracking (Nonlinear Plant)');
for kk = 1:2
    subplot(2,1,kk);
    plot(t_sp, y_sp_log(:,kk), 'b','LineWidth',1.2); hold on; grid on;
    plot(t_sp, r_sp_phys(:,kk), 'k--','LineWidth',1);
    ylabel(outNames{kk},'Interpreter','latex'); xlabel('Time [s]');
    legend('Nonlinear plant output','Reference');
end

figure('Name','Scenario A: MV Trajectories');
for i = 1:nu
    subplot(nu,1,i);
    plot(t_sp, u_sp_log(:,i),'LineWidth',1.1); hold on; grid on;
    yline(mv_bounds(i,1),'r--'); yline(mv_bounds(i,2),'r--');
    ylabel(varNames{idx_MV(i)});
end
%% Scenario B: Disturbance rejection on nonlinear plant (softer MD)
x_hat_B = x_hat;      % branch off from warmed-up state
x_nl_B  = x_nl;
u_full_B = u_full;

Tsim_d = 300;
r_d_traj = repmat(y_Normalized(1,:), Tsim_d, 1);   % hold nominal reference

% Use smaller, in-range disturbance (e.g. ±5 around nominal)
md_nom = u0(idx_MD);
d_phys_traj = repmat(md_nom, Tsim_d, 1);
d_phys_traj(50:80)  = md_nom + 5;    % small step up
d_phys_traj(120:150) = md_nom - 5;   % small step down

y_d_log = zeros(Tsim_d, ny);
u_d_log = zeros(Tsim_d, nu);

for k = 1:Tsim_d
    % Apply MD to nonlinear plant
    u_full_B(idx_MD) = d_phys_traj(k);
    md_norm_k = (u_full_B(idx_MD) - mu_u(idx_MD)) ./ sigma_u(idx_MD);

    % Measurement and correction
    y_phys = sys_output_wrapper(x_nl_B, u_full_B.', p).';
    y_norm = (y_phys - mu_y) ./ sigma_y;
    x_hat_B = x_hat_B + L * (y_norm.' - C*x_hat_B);
    mpcstate_obj.Plant = x_hat_B;

    % MPC move with measured MD
    [mv_norm, ~] = mpcmove(mpcobj, mpcstate_obj, y_norm, r_d_traj(k,:), md_norm_k);
    mv_norm = min(max(mv_norm(:), mv_bounds_norm(:,1)), mv_bounds_norm(:,2));
    mv_phys = mv_norm(:).' .* sigma_u(idx_MV) + mu_u(idx_MV);
    u_full_B(idx_MV) = mv_phys;

    % Nonlinear DAE step
    [~, x_traj] = ode15s(@(tt,xx) ode_PEMFC(tt,xx,u_full_B.'), [0 Ts], x_nl_B, options);
    x_nl_B = x_traj(end,:).';

    % Predictor step
    x_hat_B = A*x_hat_B + Bu*mv_norm(:) + Bd*md_norm_k.';
    y_d_log(k,:) = y_phys;
    u_d_log(k,:) = mv_phys;
end