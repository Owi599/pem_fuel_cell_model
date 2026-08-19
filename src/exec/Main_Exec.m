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
    ax.FontSize= 24;
    ylabel(outNames{kk},'FontSize',28,'FontWeight','bold','Interpreter','latex');
    xlabel('Time [s]','FontSize',28,'FontWeight','bold','Interpreter','latex');
    legend('Nonlinear DAE','Linear (N4SID)','FontSize',24,'FontWeight','bold','Location','best','Orientation','Horizontal');
end
data = iddata(y_Normalized, U_Normalized, Ts);
figure('Name','compare(): NRMSE fit');
compare(data, sys);
ax1= gca;
ax1.FontSize = 24;

%% Role assignment printout
fprintf('\nMV: %s\n', strjoin(varNames(idxMV), ', '));
fprintf('MD: %s\n', strjoin(varNames(idxMD), ', '));
fprintf('Fixed: %s\n', strjoin(varNames(idxUD), ', '));
nu = numel(idxMV); nd = numel(idxMD); ny = 2;
Bu = B(:, idxMV); Bd = B(:, idxMD);

%% Design MPC (default estimator — needed for linear sim tests)
mv_bounds = [ 0  534.7;
    0  127.0;
    0  235.9;
    0  139.5;
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
%% S1: Large setpoint steps toward a realistic target (80 degC, 80% RH)
T_target_phys    = 80 + 273.15;   % 353.15 K
aH2O_target_phys = 0.80;
check_target_within_envelope([T_target_phys, aH2O_target_phys], sys, mu_y, sigma_y, ...
    {'T_S','a_H2O_avg'});

Tsim_S1 = 2500;   % a step this size will likely need longer than 1500s to settle
kStep_S1 = 50;

scenarios_S1 = struct( ...
    'name', {'S1a_TempOnly','S1b_HumidityOnly','S1c_Simultaneous'}, ...
    'delta_phys', { ...
        [T_target_phys - mu_y(1), 0], ...
        [0, aH2O_target_phys - mu_y(2)], ...
        [T_target_phys - mu_y(1), aH2O_target_phys - mu_y(2)] });

resS1_all = struct();
for kS = 1:numel(scenarios_S1)
    delta_norm = scenarios_S1(kS).delta_phys ./ sigma_y;
    r_traj = zeros(Tsim_S1, 2);
    r_traj(kStep_S1:end,:) = repmat(delta_norm, Tsim_S1 - kStep_S1 + 1, 1);

    res = run_nonlinear_test(mpcobj, sys, p, options, idxMV, idxMD, ...
        mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, x0, u0, r_traj, [], Ts);

    r_phys = r_traj .* sigma_y + mu_y;
    t_vec = (0:Tsim_S1-1).' * Ts;

    plot_results(t_vec, res.y, r_phys, res.u, mv_bounds, outNames, varNames(idxMV), ...
        sprintf('%s: target [%.2f K, %.3f]', scenarios_S1(kS).name, T_target_phys, aH2O_target_phys), []);

    resS1_all.(scenarios_S1(kS).name) = res;

    % Steady-state check + spatial distribution at the end of this run
    tailFrac = round(0.9*Tsim_S1):Tsim_S1;
    y_tail_range = max(res.y(tailFrac,:)) - min(res.y(tailFrac,:));
    if any(y_tail_range > 0.02 .* abs(r_phys(end,:) - r_phys(1,:)))
        warning(['%s: outputs still moving by >2%% of the step size over the final ' ...
            '10%% of the simulation window (T_S range = %.4f K, a_H2O range = %.4f). ' ...
            'Consider extending Tsim_S1 before trusting this as a steady-state snapshot.'], ...
            scenarios_S1(kS).name, y_tail_range(1), y_tail_range(2));
    end

    x_final = res.x(end,:).';   % last logged state, transposed to column
    plot_spatial_profiles(x_final, p, ...
        sprintf('%s: Steady-State Distribution (target [%.2f K, %.3f])', ...
        scenarios_S1(kS).name, T_target_phys, aH2O_target_phys));
end
%% S2: Large setpoint steps down toward a low-operating-point target
T_target_phys_S2    = 74 + 273.15;   %  K -- adjust to your actual low target
aH2O_target_phys_S2 = 0.38;          % adjust to your actual low target

check_target_within_envelope([T_target_phys_S2, aH2O_target_phys_S2], sys, mu_y, sigma_y, ...
    {'T_S','a_H2O_avg'});

Tsim_S2 = 2500;   % match S1's longer window given the step size
kStep_S2 = 50;

scenarios_S2 = struct( ...
    'name', {'S2a_TempOnly','S2b_HumidityOnly','S2c_Simultaneous'}, ...
    'delta_phys', { ...
        [T_target_phys_S2 - mu_y(1), 0], ...
        [0, aH2O_target_phys_S2 - mu_y(2)], ...
        [T_target_phys_S2 - mu_y(1), aH2O_target_phys_S2 - mu_y(2)] });

resS2_all = struct();
for kS = 1:numel(scenarios_S2)
    delta_norm = scenarios_S2(kS).delta_phys ./ sigma_y;
    r_traj = zeros(Tsim_S2, 2);
    r_traj(kStep_S2:end,:) = repmat(delta_norm, Tsim_S2 - kStep_S2 + 1, 1);

    res = run_nonlinear_test(mpcobj, sys, p, options, idxMV, idxMD, ...
        mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, x0, u0, r_traj, [], Ts);

    r_phys = r_traj .* sigma_y + mu_y;
    t_vec = (0:Tsim_S2-1).' * Ts;

    plot_results(t_vec, res.y, r_phys, res.u, mv_bounds, outNames, varNames(idxMV), ...
        sprintf('%s: target [%.2f K, %.3f]', scenarios_S2(kS).name, T_target_phys_S2, aH2O_target_phys_S2), []);

    resS2_all.(scenarios_S2(kS).name) = res;

    % Steady-state check + spatial distribution at the end of this run
    tailFrac = round(0.9*Tsim_S2):Tsim_S2;
    y_tail_range = max(res.y(tailFrac,:)) - min(res.y(tailFrac,:));
    if any(y_tail_range > 0.02 .* abs(r_phys(end,:) - r_phys(1,:)))
        warning(['%s: outputs still moving by >2%% of the step size over the final ' ...
            '10%% of the simulation window (T_S range = %.4f K, a_H2O range = %.4f). ' ...
            'Consider extending Tsim_S2 before trusting this as a steady-state snapshot.'], ...
            scenarios_S2(kS).name, y_tail_range(1), y_tail_range(2));
    end

    x_final = res.x(end,:).';
    plot_spatial_profiles(x_final, p, ...
        sprintf('%s: Steady-State Distribution (target [%.2f K, %.3f])', ...
        scenarios_S2(kS).name, T_target_phys_S2, aH2O_target_phys_S2));
end

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
%% S4-revised: full-range current disturbance rejection
I_p05 = prctile(U(:,idxMD), 5);   % near-minimum load
I_p95 = prctile(U(:,idxMD), 95);  % near-maximum load
I_nom = u0(idxMD);                % nominal / medium load
fprintf('Current levels: min = %.3f, nominal = %.3f, max = %.3f\n', I_p05, I_nom, I_p95);

%% Dynamic transient: ramp min -> hold -> max -> hold -> nominal
Tsim_S4 = 1800;
r_S4_traj = zeros(Tsim_S4, 2);   % outputs held at nominal throughout

d_S4_phys = I_nom * ones(Tsim_S4,1);
segLen = 10;
k1 = 100; k1e = k1+segLen-1;
k2s = 600;
k3 = k2s; k3e = k3+segLen-1;
k4s = 1200;
k5 = k4s; k5e = k5+segLen-1;

d_S4_phys(k1:k1e)    = linspace(I_nom, I_p05, segLen).';
d_S4_phys(k1e+1:k2s) = I_p05;
d_S4_phys(k3:k3e)    = linspace(I_p05, I_p95, segLen).';
d_S4_phys(k3e+1:k4s) = I_p95;
d_S4_phys(k5:k5e)    = linspace(I_p95, I_nom, segLen).';
d_S4_phys(k5e+1:end) = I_nom;

resS4 = run_nonlinear_test(mpcobj, sys, p, options, idxMV, idxMD, ...
    mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, x0, u0, r_S4_traj, d_S4_phys, Ts);

r_S4_phys = r_S4_traj .* sigma_y + mu_y;
t_S4 = (0:Tsim_S4-1).' * Ts;
plot_results(t_S4, resS4.y, r_S4_phys, resS4.u, mv_bounds, outNames, varNames(idxMV), ...
    'S4: Full-Range Current Disturbance Rejection', []);

figure('Name','S4: Current profile');
plot(t_S4, d_S4_phys, 'LineWidth',1.5); grid on;
xlabel('Time [s]'); ylabel(varNames{idxMD});

%% Steady-state spatial distributions at the three requested load points
current_levels = [I_p05, I_nom, I_p95];
labels3 = {'Min load','Nominal load','Max load'};
Tsim_hold = 400;

resSteady = cell(1,3);
for kL = 1:3
    d_hold = current_levels(kL) * ones(Tsim_hold,1);
    r_hold = zeros(Tsim_hold, 2);
    resSteady{kL} = run_nonlinear_test(mpcobj, sys, p, options, idxMV, idxMD, ...
        mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, x0, u0, r_hold, d_hold, Ts);

    plot_spatial_profiles(resSteady{kL}.x(end,:).', p, ...   % note: row k, transposed to column
        sprintf('%s (I = %.3f): Steady-State Distribution', labels3{kL}, current_levels(kL)));
end
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
%%
function check_target_within_envelope(target_phys, sys, mu_y, sigma_y, labels)
y_phys_ID = sys.UserData.y_n .* sigma_y + mu_y;   % recovered physical y from ID data
for k = 1:numel(target_phys)
    lo = prctile(y_phys_ID(:,k), 5);
    hi = prctile(y_phys_ID(:,k), 95);
    if target_phys(k) < lo || target_phys(k) > hi
        warning(['%s target (%.4f) is OUTSIDE the 5-95%% identification range ' ...
            '[%.4f, %.4f]. The linear model is extrapolating here.'], ...
            labels{k}, target_phys(k), lo, hi);
    else
        fprintf('%s target (%.4f) is within ID range [%.4f, %.4f].\n', ...
            labels{k}, target_phys(k), lo, hi);
    end
end
end
%%
function plot_spatial_profiles(x_final, p, figTitle)
% x_final: full DAE state vector at (approximate) steady state.
xStruct = p.state2struct(x_final);
z = (1:p.N) * p.delta_z;

aH2O_a_z = (p.R * xStruct.T_a .* xStruct.c_H2O_a) ./ p.p_sat(xStruct.T_a);
aH2O_c_z = (p.R * xStruct.T_c .* xStruct.c_H2O_c) ./ p.p_sat(xStruct.T_c);

figure('Name', figTitle);
subplot(2,1,1);
plot(z, xStruct.T_s, 'LineWidth', 1.5); grid on;
ylabel('$T^S [K]$','Interpreter','latex','FontSize',24,'FontWeight','bold'); xlabel('Channel position $z$ [m]','Interpreter','latex','FontSize',24,'FontWeight','bold');
title([figTitle ' — Temperature distribution']);
ax1 = gca;
ax1.FontSize = 18;
subplot(2,1,2);
plot(z, aH2O_a_z, 'b', 'LineWidth', 1.5); hold on;
plot(z, aH2O_c_z, 'r', 'LineWidth', 1.5); grid on;
yline(1.0, 'k--', 'Saturation (a_{H2O} = 1)');
legend('Anode','Cathode','Location','best','Interpreter','latex','FontSize',20,'FontWeight','bold','Orientation','Horizontal');
ylabel('$a_{H2O}$ [-]','Interpreter','latex','FontSize',24,'FontWeight','bold'); xlabel('Channel position $z$ [m]','Interpreter','latex','FontSize',24,'FontWeight','bold');
title([figTitle ' — Water activity distribution'],'Interpreter','latex','FontSize',24,'FontWeight','bold');
ax2 = gca;
ax2.FontSize = 18;
if any(aH2O_a_z > 0.98) || any(aH2O_c_z > 0.98)
    warning(['%s: water activity approaches/exceeds saturation somewhere in the ' ...
        'channel. This vapor-only model does not capture liquid water formation ' ...
        'reliably — treat this region with caution.'], figTitle);
end
end
%% S7: Plant parameter perturbation (model-plant mismatch)
% Purpose: sys/mpcobj were identified around ONE operating condition of
% the nominal plant p. Here we perturb the *nonlinear plant* parameters
% while keeping the controller (mpcobj, built from the nominal sys)
% fixed, to see how much mismatch the closed loop tolerates.
%
% p exposes exactly this kind of knob already: p.k_kappa, p.k_D_w,
% p.k_lambda, p.k_t are dimensionless multiplicative gains on membrane
% conductivity, water diffusivity, water-uptake isotherm, and
% electro-osmotic drag/transport number respectively (see the
% "additional gains" block in mod_param_PEMFC). Nominally all = 1, so a
% fractional perturbation here is exactly (1 + pertLevel).
%
% Two more physically-grounded scalar knobs worth including:
%   p.alpha_1 / p.alpha_2  - solid<->gas / coolant<->solid heat transfer
%                             coefficients (directly drives T^S dynamics)
%   p.i_0_ref_c             - cathode reference exchange current density
%                             (ORR kinetics; strong effect on U_cell and
%                             indirectly on the thermal balance)

pertLevels = [-0.20, -0.10, 0, +0.10, +0.20];   % fractional perturbation
paramsToPerturb = {'k_kappa', 'k_D_w', 'alpha_1', 'i_0_ref_c'};

Tsim_S7 = 800;
r_S7_traj = zeros(Tsim_S7, 2);   % nominal regulation target

resS7 = struct();
for iParam = 1:numel(paramsToPerturb)
    fname = paramsToPerturb{iParam};
    if ~isfield(p, fname)
        warning('S7: parameter "%s" not found on p; skipping.', fname);
        continue;
    end
    nominalVal = p.(fname);

    for iLevel = 1:numel(pertLevels)
        p_pert = p;                                  % shallow copy
        p_pert.(fname) = nominalVal * (1 + pertLevels(iLevel));

        res = run_nonlinear_test(mpcobj, sys, p_pert, options, idxMV, idxMD, ...
            mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, x0, u0, ...
            r_S7_traj, [], Ts);

        tag = sprintf('%s_%+d pct', fname, round(100*pertLevels(iLevel)));
        resS7.(matlab.lang.makeValidName(tag)) = res;

        % Quantify degradation: RMSE of the tracking error over the run
        % (res.y is in normalized deviation coords; multiply by sigma_y
        % for physical units, e.g. rmse_T_K = rmse_T .* sigma_y(1))
        rmse_T   = sqrt(mean((res.y(:,1)).^2));
        rmse_aH2 = sqrt(mean((res.y(:,2)).^2));
        fprintf('S7 [%s]: RMSE_T=%.4f, RMSE_aH2O=%.5f\n', tag, rmse_T, rmse_aH2);
    end
end
% Suggested follow-up: plot RMSE vs. pertLevels per parameter (a
% tornado/sensitivity chart) to see which of {k_kappa, k_D_w, alpha_1,
% i_0_ref_c} the controller is most fragile to. i_0_ref_c is a good
% candidate to show the biggest sensitivity since it enters the
% electrochemical source terms directly, not just a transport gain.

%% S8: Aged / degraded stack test
% Purpose: single, larger, more "real" mismatch case: an aged stack
% (reduced performance) controlled by a model identified on a fresh
% stack. Degradation mechanisms modeled here (all physically motivated,
% magnitudes are illustrative -- tune to your target aging state):
%   - catalyst activity loss (ORR):      i_0_ref_c  -30%
%   - membrane conductivity loss
%     (drying / chemical degradation):   k_kappa    -15%
%   - cathode flow-channel restriction
%     (GDL degradation / flooding):      K_c        -20%

p_aged = p;
p_aged.i_0_ref_c = p.i_0_ref_c * 0.70;
p_aged.k_kappa    = p.k_kappa    * 0.85;
p_aged.K_c        = p.K_c        * 0.80;

Tsim_S8 = 1000;
r_S8_traj = zeros(Tsim_S8, 2);

resS8 = run_nonlinear_test(mpcobj, sys, p_aged, options, idxMV, idxMD, ...
    mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, x0, u0, r_S8_traj, [], Ts);

r_S8_phys = r_S8_traj .* sigma_y + mu_y;
t_S8 = (0:Tsim_S8-1).' * Ts;
plot_results(t_S8, resS8.y, r_S8_phys, resS8.u, mv_bounds, outNames, varNames(idxMV), ...
    'S8: Aged Stack Regulation (Model-Plant Mismatch)', []);

%% S9: Sensor noise robustness
% Purpose: inject measurement noise before it reaches the estimator and
% check tracking degradation + control-effort chatter.
%
% NOTE: this requires access to the measurement injection point inside
% run_nonlinear_test. If run_nonlinear_test does not currently support a
% noise input, the minimal change is to add an optional argument
% 'measNoiseStd' (a 1x2 vector, one std per output, in PHYSICAL units)
% and add it to y_meas = y_true + measNoiseStd .* randn(size(y_true))
% right before the measurement is fed to the estimator update step.

rng(42);  % reproducibility
measNoiseStd_phys = [0.05, 0.002];   % [K, water activity] 1-sigma sensor noise -- EDIT to your sensor specs

Tsim_S9 = 600;
r_S9_traj = repmat(y_Normalized(1,:), Tsim_S9, 1);

% If run_nonlinear_test supports a noise argument:
% resS9 = run_nonlinear_test(mpcobj, sys, p, options, idxMV, idxMD, ...
%     mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, x0, u0, ...
%     r_S9_traj, [], Ts, [], [], measNoiseStd_phys);   % <-- extend signature

% Fallback if signature can't be changed: run several seeds and compare
% output variance/control effort against the noiseless S0/base run to
% at least gauge sensitivity qualitatively via repeated identical runs
% (won't show noise rejection directly without a real injection point).
nSeeds = 5;
u_effort = zeros(nSeeds,1);
for iSeed = 1:nSeeds
    rng(iSeed);
    resTmp = run_nonlinear_test(mpcobj, sys, p, options, idxMV, idxMD, ...
        mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, x0, u0, ...
        r_S9_traj, [], Ts);
    du = diff(resTmp.u, 1, 1);
    u_effort(iSeed) = sum(du(:).^2);   % control-move "energy" as a chatter proxy
end
fprintf('S9 (no-noise-injection fallback): control effort mean=%.4g, std=%.4g\n', ...
    mean(u_effort), std(u_effort));
% Once measurement noise injection exists in run_nonlinear_test, replace
% this loop with the commented call above and plot y-tracking + du(t).

%% S10: Sensor bias / drift
% Purpose: slowly-drifting bias on one measured output; check how long
% the controller "happily" regulates to the wrong physical value.
% Same prerequisite as S9: run_nonlinear_test needs a hook to add a bias
% to the measured (not true) output before the estimator sees it.

Tsim_S10 = 900;
biasRampStart = 100;
biasFinal_phys = [0.3, 0];    % +0.3 K drift on T^S only -- EDIT
biasProfile = zeros(Tsim_S10, 2);
biasProfile(biasRampStart:end, 1) = linspace(0, biasFinal_phys(1), ...
    Tsim_S10 - biasRampStart + 1);

r_S10_traj = zeros(Tsim_S10, 2);

% resS10 = run_nonlinear_test(mpcobj, sys, p, options, idxMV, idxMD, ...
%     mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, x0, u0, ...
%     r_S10_traj, [], Ts, [], [], [], biasProfile);   % <-- extend signature
%
% After running, compare resS10.y (what the controller believed) against
% the TRUE plant output logged separately (you may need
% run_nonlinear_test to also return the unbiased/true y for this
% comparison to be meaningful) to quantify the steady-state physical
% offset the controller settles into.

%% S11: Stuck actuator (MV fault)
% Purpose: freeze one MV at a fixed value partway through a run and see
% whether the remaining MVs can compensate, and how much the coupled
% output degrades.

Tsim_S11 = 800;
kFault_S11 = 200;
stuckMVname = 'n_dot_H2O_a_in';   % EDIT to the MV you want to fail
iStuckMV = find(strcmp(varNames(idxMV), stuckMVname));
assert(isscalar(iStuckMV), 'S11: could not locate stuck MV "%s".', stuckMVname);

mpcS11 = mpcobj;
u0_norm_mv = (u0(idxMV) - mu_u(idxMV)) ./ sigma_u(idxMV);
stuckVal_norm = u0_norm_mv(iStuckMV);   % freeze at its nominal value

% Emulate "stuck" by collapsing that MV's Min/Max to the frozen value
% from kFault_S11 onward. If run_nonlinear_test runs the sim in one
% shot, split the simulation into two segments: pre-fault (normal) and
% post-fault (with the tightened mpc object), chaining plant/estimator
% state between them.

r_S11_traj = repmat([0.5, 0] ./ sigma_y, Tsim_S11, 1);   % modest T step, requested throughout

% --- Segment 1: pre-fault, normal controller ---
resS11_pre = run_nonlinear_test(mpcobj, sys, p, options, idxMV, idxMD, ...
    mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, x0, u0, ...
    r_S11_traj(1:kFault_S11,:), [], Ts);

% --- Segment 2: post-fault, MV frozen ---
mpcS11.MV(iStuckMV).Min = stuckVal_norm;
mpcS11.MV(iStuckMV).Max = stuckVal_norm;
mpcS11.MV(iStuckMV).RateMin = 0;
mpcS11.MV(iStuckMV).RateMax = 0;

% Chain initial conditions from the end of segment 1.
% TODO: confirm run_nonlinear_test exposes the final plant state (x)
% and controller/estimator state needed to warm-start segment 2 -- if
% it currently only accepts x0 (plant DAE state) as in the main script,
% pass resS11_pre.x(end,:).' here; if it also needs an MPC/estimator
% state carry-over, extend the function to accept/return one.
x0_seg2 = resS11_pre.x(end,:).';

resS11_post = run_nonlinear_test(mpcS11, sys, p, options, idxMV, idxMD, ...
    mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, x0_seg2, u0, ...
    r_S11_traj(kFault_S11+1:end,:), [], Ts);

% Stitch results for plotting.
y_S11 = [resS11_pre.y; resS11_post.y];
u_S11 = [resS11_pre.u; resS11_post.u];
r_S11_phys = r_S11_traj .* sigma_y + mu_y;
t_S11 = (0:Tsim_S11-1).' * Ts;

plot_results(t_S11, y_S11, r_S11_phys, u_S11, mv_bounds, outNames, varNames(idxMV), ...
    sprintf('S11: Stuck Actuator Fault (%s frozen at t=%d s)', stuckMVname, kFault_S11), []);

%% S12: Step (non-ramped) current disturbance -- worst-case transient
% Purpose: your existing S4 ramps the current over 10 samples; here we
% apply an instantaneous step to probe worst-case overshoot / MV
% saturation on the same min/nominal/max load levels.

I_p05_S12 = prctile(U(:,idxMD), 5);
I_p95_S12 = prctile(U(:,idxMD), 95);
I_nom_S12 = u0(idxMD);

Tsim_S12 = 1200;
r_S12_traj = zeros(Tsim_S12, 2);   % hold outputs at nominal

d_S12_phys = I_nom_S12 * ones(Tsim_S12,1);
kStepUp   = 100;
kStepDown = 600;
d_S12_phys(kStepUp:kStepDown-1) = I_p95_S12;    % instantaneous step to max load
d_S12_phys(kStepDown:end)       = I_p05_S12;    % instantaneous step to min load

resS12 = run_nonlinear_test(mpcobj, sys, p, options, idxMV, idxMD, ...
    mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, x0, u0, ...
    r_S12_traj, d_S12_phys, Ts);

r_S12_phys = r_S12_traj .* sigma_y + mu_y;
t_S12 = (0:Tsim_S12-1).' * Ts;

plot_results(t_S12, resS12.y, r_S12_phys, resS12.u, mv_bounds, outNames, varNames(idxMV), ...
    'S12: Step (Non-Ramped) Current Disturbance -- Worst-Case Transient', []);

figure('Name','S12: Current profile (step)');
plot(t_S12, d_S12_phys, 'LineWidth',1.5); grid on;
xlabel('Time [s]'); ylabel(varNames{idxMD});

% Report peak MV excursion vs. bounds -- flags saturation immediately.
for iMV = 1:numel(idxMV)
    uPhys_iMV = resS12.u(:,iMV);
    fprintf('S12 MV "%s": min=%.3f, max=%.3f, bounds=[%.3f, %.3f]\n', ...
        varNames{idxMV(iMV)}, min(uPhys_iMV), max(uPhys_iMV), ...
        mv_bounds(iMV,1), mv_bounds(iMV,2));
end
