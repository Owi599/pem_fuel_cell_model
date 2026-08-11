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

outNames = {'$T^S(z = 123.1 mm) [K]$','$\overline{a}_{H2O} [-]$'};
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
delta_r_phys = [2, 0.02];

r1_phys = mu_y + delta_r_phys;
r1 = (r1_phys - mu_y) ./ sigma_y;

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


%% Use MPC Toolbox built-in Kalman estimator for nonlinear tests

[Lmpc, Mmpc, Aest, Cmest, Buest] = getEstimator(mpcobj);

options.Mass = p.M();
options.RelTol = 1e-6;
options.AbsTol = 1e-6;
options.MStateDependence = 'none';


%% Closed-loop nominal test (base case)
Tsim_nl = 500;
r_nl_traj = repmat(y_Normalized(1,:), Tsim_nl, 1);
resBase = run_nonlinear_test(mpcobj, sys, p, options,  idxMV, idxMD, ...
    mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, x0, u0 , ...
    r_nl_traj, [], Ts);


%% Scenario S0: Nominal regulation on nonlinear plant
Tsim_S0 = 600;

% Zero deviations in normalized coordinates:
% desired physical outputs remain equal to the nominal operating point.
r_S0_traj = zeros(Tsim_S0, 2);

resS0 = run_nonlinear_test( ...
    mpcobj, sys, p, options,  idxMV, idxMD, ...
    mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, ...
    x0, u0,  r_S0_traj, [], Ts);

% Convert the normalized reference to physical units.
r_S0_phys = r_S0_traj .* sigma_y + mu_y;

% Physical time vector.
t_S0 = (0:Tsim_S0-1).' * Ts;
ytickformat('%.2f');

% Half-width of the displayed physical y-axis ranges:
% [temperature in K, water activity]
yHalfRange_S0 = [0.05, 0.005];

plot_results( ...
    t_S0, resS0.y, r_S0_phys, resS0.u, mv_bounds, ...
    outNames, varNames(idxMV), ...
    'Scenario S0: Nominal Regulation (Nonlinear Plant)', ...
    yHalfRange_S0);
%% S1: Positive simultaneous setpoint step
Tsim_S1 = 1500;
kStep_S1 = 50;

% Physical setpoint deviations from the nominal operating point.
deltaS1_phys = [1.0, 0.005];   % [K, water activity]

% Normalized reference trajectory used by the MPC.
r_S1_traj = zeros(Tsim_S1, 2);
r_S1_traj(kStep_S1:end,:) = repmat( ...
    deltaS1_phys ./ sigma_y, ...
    Tsim_S1 - kStep_S1 + 1, 1);

resS1 = run_nonlinear_test( ...
    mpcobj, sys, p, options, idxMV, idxMD, ...
    mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, ...
    x0, u0, r_S1_traj, [], Ts);

% Convert reference from normalized deviations to physical units.
r_S1_phys = r_S1_traj .* sigma_y + mu_y;
t_S1 = (0:Tsim_S1-1).' * Ts;

plot_results( ...
    t_S1, resS1.y, r_S1_phys, resS1.u, mv_bounds, ...
    outNames, varNames(idxMV), ...
    'S1: Positive Simultaneous Setpoint Step',[]);
%% S2: Negative simultaneous setpoint step
Tsim_S2 = 1500;
kStep_S2 = 50;

% Equal-magnitude negative reference change.
deltaS2_phys = [-1.0, -0.005];  % [K, water activity]

r_S2_traj = zeros(Tsim_S2, 2);
r_S2_traj(kStep_S2:end,:) = repmat( ...
    deltaS2_phys ./ sigma_y, ...
    Tsim_S2 - kStep_S2 + 1, 1);

resS2 = run_nonlinear_test( ...
    mpcobj, sys, p, options,  idxMV, idxMD, ...
    mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, ...
    x0, u0,  r_S2_traj, [], Ts);

r_S2_phys = r_S2_traj .* sigma_y + mu_y;
t_S2 = (0:Tsim_S2-1).' * Ts;
plot_results( ...
    t_S2, resS2.y, r_S2_phys, resS2.u, mv_bounds, ...
    outNames, varNames(idxMV), ...
    'S2: Negative Simultaneous Setpoint Step',[]);
%% S3: Sequential MIMO setpoint tracking and coupling test
Tsim_S3 = 1200;

% Physical target deviations from the nominal operating point.
deltaT_phys  = 1.0;     % K
deltaAH_phys = 0.005;   % water activity [-]

deltaT_norm  = deltaT_phys  / sigma_y(1);
deltaAH_norm = deltaAH_phys / sigma_y(2);

% Reference in normalized deviation coordinates.
r_S3_traj = zeros(Tsim_S3, 2);

% 0–49 s: nominal regulation [0, 0].
%
% 50–249 s: request temperature only, keep humidity nominal.
r_S3_traj(50:249,1) = deltaT_norm;

% 250–499 s: retain temperature target and add humidity target.
r_S3_traj(250:499,:) = repmat( ...
    [deltaT_norm, deltaAH_norm], ...
    250, 1);

% 500–699 s: humidity back to nominal; retain temperature target.
r_S3_traj(500:699,1) = deltaT_norm;

% 700–900 s: both outputs return to nominal [0, 0].

resS3 = run_nonlinear_test( ...
    mpcobj, sys, p, options,  idxMV, idxMD, ...
    mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, ...
    x0, u0,  r_S3_traj, [], Ts);

% Physical-unit reference and time vector for plots.
r_S3_phys = r_S3_traj .* sigma_y + mu_y;
t_S3 = (0:Tsim_S3-1).' * Ts;

plot_results( ...
    t_S3, resS3.y, r_S3_phys, resS3.u, mv_bounds, ...
    outNames, varNames(idxMV), ...
    'S3: Sequential Setpoint Tracking and MIMO Interaction',[]);
%% S4: Measured-disturbance rejection on nonlinear PEMFC plant
Tsim_S4 = 700;

% Keep both controlled outputs at their nominal operating-point values.
r_S4_traj = zeros(Tsim_S4, 2);

% Nominal physical value of the measured disturbance channel.
md_nom = u0(idxMD);

% Use an amplitude from the data range to stay in a realistic local region.
md_p05 = prctile(U(:,idxMD), 5);
md_p95 = prctile(U(:,idxMD), 95);
md_span = md_p95 - md_p05;

% Start conservatively: 5% of the 5th-to-95th percentile data range.
dAmp = 0.05 * md_span;

% Do not exceed the available trajectory envelope.
md_up = min(md_nom + dAmp, md_p95);

% Physical measured-disturbance trajectory.
d_S4_phys = md_nom * ones(Tsim_S4, 1);

% Smooth ramp avoids an unrealistically discontinuous DAE input.
kStart = 150;
rampLength = 10;
kRampEnd = kStart + rampLength - 1;
kHoldEnd = 400;
kReturnEnd = kHoldEnd + rampLength;

d_S4_phys(kStart:kRampEnd) = linspace(md_nom, md_up, rampLength).';
d_S4_phys(kRampEnd+1:kHoldEnd) = md_up;
d_S4_phys(kHoldEnd+1:kReturnEnd) = ...
    linspace(md_up, md_nom, rampLength).';

% The remaining samples are already nominal.
resS4 = run_nonlinear_test( ...
    mpcobj, sys, p, options,  idxMV, idxMD, ...
    mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, ...
    x0, u0,  r_S4_traj, d_S4_phys, Ts);

% Convert the zero normalized references to physical output references.
r_S4_phys = r_S4_traj .* sigma_y + mu_y;
t_S4 = (0:Tsim_S4-1).' * Ts;

plot_results( ...
    t_S4, resS4.y, r_S4_phys, resS4.u, mv_bounds, ...
    outNames, varNames(idxMV), ...
    'S4: Measured-Disturbance Rejection (Nonlinear Plant)',[]);
%% S5: Unmeasured-disturbance rejection on nonlinear PEMFC plant
Tsim_S5 = 900;

% Maintain both outputs at the nominal physical operating point.
r_S5_traj = zeros(Tsim_S5, 2);

% Select one physically meaningful unmeasured-input channel.
jUD = idxUD(1);     % Inspect varNames(jUD) and choose deliberately.

ud_nom = u0(jUD);

% Use its identification-data envelope to choose a local, realistic step.
ud_p05 = prctile(U(:,jUD), 5);
ud_p95 = prctile(U(:,jUD), 95);
ud_span = ud_p95 - ud_p05;

% Begin with +5% of its central data range.
ud_amp = 0.05 * ud_span;
ud_step = min(ud_nom + ud_amp, ud_p95);

% Plant-only, physical disturbance profile:
% ramp 150–159 s; hold until 400 s; return to nominal by 410 s.
ud_S5_phys = ud_nom * ones(Tsim_S5, 1);

kStart = 150;
rampLength = 10;
kRampEnd = kStart + rampLength - 1;
kHoldEnd = 400;
kReturnEnd = kHoldEnd + rampLength;

ud_S5_phys(kStart:kRampEnd) = linspace(ud_nom, ud_step, rampLength).';
ud_S5_phys(kRampEnd+1:kHoldEnd) = ud_step;
ud_S5_phys(kHoldEnd+1:kReturnEnd) = ...
    linspace(ud_step, ud_nom, rampLength).';

% Important: idxUD and ud_S5_phys are supplied only to the nonlinear plant.
resS5 = run_nonlinear_test( ...
    mpcobj, sys, p, options, ...
    idxMV, idxMD, mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, ...
    x0, u0, r_S5_traj, [], Ts, ...
    jUD, ud_S5_phys);

r_S5_phys = r_S5_traj .* sigma_y + mu_y;
t_S5 = (0:Tsim_S5-1).' * Ts;

plot_results( ...
    t_S5, resS5.y, r_S5_phys, resS5.u, mv_bounds, ...
    outNames, varNames(idxMV), ...
    sprintf('S5: Unmeasured Disturbance Rejection (%s)', varNames{jUD}),[]);
%% S6: Constraint-handling test with restricted humidification authority
Tsim_S6 = 700;
kStep_S6 = 50;

% Identify the LOCAL MV positions in mv_bounds / mpcobj.MV.
% Based on idxMV = [1 2 3 4 8]:
% n_dot_H2O_a_in is local MV 2
% n_dot_H2O_c_in is local MV 4
iH2Oa = find(strcmp(varNames(idxMV), 'n_dot_H2O_a_in'));
iH2Oc = find(strcmp(varNames(idxMV), 'n_dot_H2O_c_in'));

assert(isscalar(iH2Oa) && isscalar(iH2Oc), ...
    'Could not uniquely locate the two humidification MVs.');

% Copy the normal controller; retain all its tuning.
mpcS6 = mpcobj;

% Restrict only upward humidification authority:
% allow each humidification flow to rise only 10% of its original
% available nominal-to-upper-bound margin.
u0_norm_mv = (u0(idxMV) - mu_u(idxMV)) ./ sigma_u(idxMV);

alpha = 0.10;   % 10% available positive headroom

for iMV = [iH2Oa, iH2Oc]

    oldMax = mpcS6.MV(iMV).Max;
    nominal = u0_norm_mv(iMV);

    % New maximum remains above nominal, but is intentionally restrictive.
    newMax = nominal + alpha * (oldMax - nominal);

    mpcS6.MV(iMV).Max = newMax;

    % Keep the rate bound consistent with the newly available range.
    % This is optional but avoids a large first control move.
    mpcS6.MV(iMV).RateMax = min( ...
        mpcS6.MV(iMV).RateMax, ...
        0.10 * (newMax - nominal));
end

% Keep physical limits hard: a physical actuator must not exceed its limit.
mpcS6.MV(iH2Oa).MaxECR = 0;
mpcS6.MV(iH2Oc).MaxECR = 0;

review(mpcS6);
% Reference trajectory in normalized deviation coordinates.
r_S6_traj = zeros(Tsim_S6, 2);

deltaS6_phys = [0.0, 0.010];  % [K, water activity]
deltaS6_norm = deltaS6_phys ./ sigma_y;

r_S6_traj(kStep_S6:end,:) = repmat( ...
    deltaS6_norm, Tsim_S6 - kStep_S6 + 1, 1);

resS6 = run_nonlinear_test( ...
    mpcS6, sys, p, options,  idxMV, idxMD, ...
    mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, ...
    x0, u0,  r_S6_traj, [], Ts);

r_S6_phys = r_S6_traj .* sigma_y + mu_y;
t_S6 = (0:Tsim_S6-1).' * Ts;

% Convert the tightened normalized bounds back to physical units
% so the plotting function shows the actual S6 constraints.
mv_bounds_S6 = mv_bounds;

for iMV = 1:numel(idxMV)
    mv_bounds_S6(iMV,1) = mu_u(idxMV(iMV)) + ...
        sigma_u(idxMV(iMV)) * mpcS6.MV(iMV).Min;

    mv_bounds_S6(iMV,2) = mu_u(idxMV(iMV)) + ...
        sigma_u(idxMV(iMV)) * mpcS6.MV(iMV).Max;
end

plot_results( ...
    t_S6, resS6.y, r_S6_phys, resS6.u, mv_bounds_S6, ...
    outNames, varNames(idxMV), ...
    'S6: Constrained Humidity Tracking (Nonlinear Plant)',[]);