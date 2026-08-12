%% Testing file
close all; clear; clc;
ensure_ode_pemfc_fresh();
%% Load dataset
inputFileName = 'dataset_100_251117_v06';
initial_input = [2.1, 560, 30, 70, 2.3, 164.73, 30, 70, 70, 0, 426, 2.7, 2.4];
Ts = 1;
%% Generate normalized identification data via the existing pipeline
[sys12, x0, x0_est, u0, mu_u, sigma_u, mu_y, sigma_y, U_n, y_n, t, idxMV, idxMD, idxUD, p, U] = ...
    linearize_pemfc_N4SID(inputFileName, initial_input);
assert(abs(mean(diff(t)) - Ts) < 1e-10, 'Ts mismatch between exec script and identified data.');
%% Package as iddata for the sweep
zFull = iddata(y_n, U_n, Ts);
zFull.OutputName = {'T_S','a_H2O_avg'};
zFull.InputName  = p.inputs.variableNames;
%% Run the order/horizon sweep
[results, sysBest, best, diagData] = select_n4sid_model(zFull, ...
    'SplitFraction', 0.70, ...
    'OrderCand', 4:2:16, ...
    'OutputWeight', diag([0.1, 1]), ...
    'FitWeights', [0.5, 0.5], ...   % adjust if you want to weight y1/y2 unevenly
    'WriteCSV', true, ...
    'CSVFile', sprintf('%s_n4sid_sweep.csv', inputFileName), ...
    'MakeFigures', true);

%% Selected model: pull diagnostics already computed by the sweep
fprintf('\nSelected model: order %d, N4Horizon = [%d %d %d]\n', ...
    best.Order, best.r, best.sy, best.su);
fprintf('Fit y1 = %.2f %% | Fit y2 = %.2f %% | Score = %.3f\n', ...
    best.Fit_y1, best.Fit_y2, best.Score);
fprintf('Time constants: %.2f s (fastest) to %.2f s (slowest)\n', ...
    best.MinTimeConst_s, best.MaxTimeConst_s);

if best.HasNegRealPole
    warning('Selected model has a discrete pole near the negative real axis (fast, likely non-physical mode). Consider the next-best candidate instead.');
end

%% MPC diagnostic pass
[Lmpc, Mmpc, Aest, Cmest, Buest] = getEstimator(mpcobj);

% 1. Is the augmented estimator still detectable, or did the disturbance
%    model get silently reduced? review() prints this, but check pole
%    count explicitly too.
review(mpcobj);
fprintf('Augmented estimator state dimension: %d (plant order: %d)\n', ...
    size(Aest,1), size(mpcobj.Model.Plant.A,1));
% If Aest is only plant-order (no extra states beyond nx), the output
% disturbance model got dropped -- that alone would explain offset.

% 2. MV saturation check on your existing scenario results
%    (repeat for whichever scenario shows the residual error, e.g. S1/S3)
figure('Name','MV Saturation Check');
for i = 1:numel(idxMV)
    subplot(numel(idxMV),1,i);
    plot(t1, u1(:,i)); hold on; grid on;
    yline(mv_bounds_norm(i,1),'r--'); yline(mv_bounds_norm(i,2),'r--');
    title(sprintf('MV %d: %s', i, varNames{idxMV(i)}));
end

% 3. Terminal tracking error, numeric not just visual
e_final = r1(end,:) - y1(end,:);
fprintf('Final tracking error: T_S = %.4f, a_H2O = %.4f (normalized)\n', e_final);
% Also check the trend of the last ~20% of the trajectory:
tailFrac = round(0.8*numel(t1)):numel(t1);
e_tail = r1 - y1(tailFrac,:);
fprintf('Mean |error| over final 20%% of window: %.4f, %.4f\n', mean(abs(e_tail)));
%% Corrected nonlinear diagnostics (physical units throughout)
r_check = r_S1_phys;      % physical reference, not r_S1_traj
y_check = resS1.y;        % already physical
u_check = resS1.u;        % already physical
t_check = t_S1;

e_final = r_check(end,:) - y_check(end,:);
fprintf('Nonlinear plant final tracking error: T_S = %.4f K, a_H2O = %.4f\n', e_final);

tailFrac = round(0.8*numel(t_check)):numel(t_check);
e_tail = r_check(tailFrac,:) - y_check(tailFrac,:);
fprintf('Nonlinear plant mean |error| over final 20%%: %.4f K, %.4f\n', mean(abs(e_tail)));

figure('Name','Nonlinear MV Saturation Check (S1)');
for i = 1:numel(idxMV)
    subplot(numel(idxMV),1,i);
    plot(t_check, u_check(:,i)); hold on; grid on;
    yline(mv_bounds(i,1),'r--'); yline(mv_bounds(i,2),'r--');   % physical bounds
    title(sprintf('MV %d: %s', i, varNames{idxMV(i)}));
end