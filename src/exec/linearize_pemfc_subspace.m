%% PEMFC: N4SID linear model 
close all; clear; clc;

%% Setup
ensure_ode_pemfc_fresh(); % only run if something was changed in mod_param_PEMFC.m

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

%normalize_u = @(U) (U - mu_u) ./ sigma_u;
%normalize_y = @(Y) (Y - mu_y) ./ sigma_y;
%U_Normalized = normalize_u(U);
%y_Normalized = normalize_y(y);
fprintf('\nPer-channel training stats (mean, std):\n');
fprintf('  y1 (T_S(10))    : mean=%.3g, std=%.3g\n', mu_y(1), sigma_y(1));
fprintf('  y2 (a_H2O_avg)  : mean=%.3g, std=%.3g\n', mu_y(2), sigma_y(2));

%% Identify linear model (N4SID, raw physical units) 
% The issue is mainly in this section, I can't get a fitting higher than
% 52% for y2 (a_H2O_AVG)

Ts = 1;
nx = 10:30;
horizons = [5 15 15];

data = iddata(y,U,Ts);

%OPT = ssestOptions( 'Focus','simulation', ...
%  'EnforceStability', false);
%OPT.OutputWeight =diag([0.1,1]);
opt = n4sidOptions('N4Weight','CVA','Focus','simulation', ...
    'N4Horizon','auto', 'EnforceStability',false, ...
    'OutputWeight', diag(1./sigma_y.^2));
[sys, x0_est] = n4sid(U, y, nx, 'Ts', Ts, opt);
A = sys.A; B = sys.B; C = sys.C; D = sys.D;
nx_id = size(A,1);
%sys_refined = ssest(data,sys,OPT);  % refining the system

fprintf('Stable:       %d\n', all(abs(eig(A)) < 1));
fprintf('Controllable: %d\n', rank(ctrb(A,B)) == nx_id);
fprintf('Observable:   %d\n', rank(obsv(A,C)) == nx_id);
fprintf('Model order:  %d\n', nx_id);

% Linear vs nonlinear open-loop comparison

y_lin = lsim(sys, U, t,x0_est);
outNames = {'T_S(10) [K]','a_H2O_{avg} [-]'};
figure('Name','Linear (N4SID) vs Nonlinear (DAE) open-loop response');
for kk = 1:2
    subplot(2,1,kk);
    plot(t, y(:,kk), 'g', 'LineWidth', 1.2); hold on;
    plot(t, y_lin(:,kk), 'r--', 'LineWidth', 1.2);
    ylabel(outNames{kk}); xlabel('Time [s]');
    legend('Nonlinear DAE','Linear (N4SID)');
end
figure('Name','compare(): NRMSE fit, linear model vs nonlinear DAE data');
compare(iddata(y, U, Ts), sys);

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

mv_bounds = [ 81.0  534.7;
               1.03  127.0;
              26.9  235.9;
               0.92  139.5;
             305.0  352.7];
range_u = diff(mv_bounds, 1, 2);

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
mpcobj.Weights.OutputVariables          = [1 1];
mpcobj.Weights.ManipulatedVariables     = 0.01 * ones(1, numel(idx_MV));
mpcobj.Weights.ManipulatedVariablesRate = 0.1  * ones(1, numel(idx_MV));

% No manual ScaleFactor tuning at this stage
% (let MPC auto-scale internally)
for i = 1:numel(idx_MV)
    mpcobj.MV(i).Min = mv_bounds(i,1);
    mpcobj.MV(i).Max = mv_bounds(i,2);
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
p_horizon = 30;   % Prediction horizon
m_horizon = 5;    % Control horizon
mpcobj = mpc(plant, Ts, p_horizon, m_horizon);

mpcobj.Weights.OutputVariables = 1 ./ sigma_y;
mpcobj.Weights.ManipulatedVariables      = 0.01 * ones(1, nu);    % small usage penalty
mpcobj.Weights.ManipulatedVariablesRate  = 0.1 ./ range_u(:)';    % keep these
mpcobj.Weights.ManipulatedVariablesRate = 20 * mpcobj.Weights.ManipulatedVariablesRate;
for i = 1:nu
    mpcobj.MV(i).ScaleFactor = range_u(i);
end
mpcobj.OV(1).ScaleFactor = sigma_y(1);
mpcobj.OV(2).ScaleFactor = sigma_y(2);
for i = 1:numel(idx_MV)
    mpcobj.MV(i).Min = mv_bounds(i,1);
    mpcobj.MV(i).Max = mv_bounds(i,2);
end
review(mpcobj);             % sanity check
mpcstate_obj = mpcstate(mpcobj);

%% Set point tracking on linear system
Tsim1 = 2000;
r1 = [mu_y(1)+2*sigma_y(1), mu_y(2)+0.5*sigma_y(2)];  % step change in both refs
r1_traj = repmat(r1, Tsim1, 1);
simopt = mpcsimopt(mpcobj);
simopt.PlantInitialState = x0_est;     % nominal linear state matching mu_y
simopt.MVSignal = [];                  % leave default unless you want open-loop MV override

[y1, t1, u1] = sim(mpcobj, Tsim1, r1_traj, [], simopt);

figure('Name','Test 1: Setpoint Tracking');
for kk = 1:2
    subplot(2,1,kk);
    plot(t1, y1(:,kk), 'b','LineWidth',1.2); hold on;
    yline(r1(kk),'k--');
    ylabel(outNames{kk}); xlabel('Time [s]');
    legend('Output','Reference');
end

%% Disturbance rejection on linear system
Tsim2 = 2000;
r2_traj = repmat(mu_y, Tsim2, 1);        % hold at nominal steady output
d2_traj = zeros(Tsim2,1);
d2_traj(50:end) = mu_u(idx_MD) + 2*sigma_u(idx_MD);  % step disturbance in I_cell at t=50

[y2s, t2, u2s] = sim(mpcobj, Tsim2, r2_traj, d2_traj);

figure('Name','Test 2: Disturbance Rejection');
for kk = 1:2
    subplot(2,1,kk);
    plot(t2, y2s(:,kk), 'r','LineWidth',1.2); hold on;
    yline(mu_y(kk),'g--');
    ylabel(outNames{kk}); xlabel('Time [s]');
    legend('Output','Nominal Reference');
end

