%% sweep_mpc_settings.m
% ------------------------------------------------------------------
% Sweeps MPC tuning parameters (prediction horizon, control horizon
% blocking, MV rate weights, OV weights) on the linear PEMFC model,
% then tests each variant in closed loop on the NONLINEAR plant and
% scores it on tracking error, overshoot, settling time, actuator
% effort, and saturation.
%
% ASSUMPTIONS / PREREQS
% ----------------------------------------------------------------
% Run this AFTER the first part of your Main_Exec.m, i.e. after:
%   [sys, x0, x0_est, u0, mu_u, sigma_u, mu_y, sigma_y, ...
%       U_Normalized, y_Normalized, t, idxMV, idxMD, idxUD, p, U] = ...
%       linearize_pemfc_N4SID(...)
% so the following variables already exist in the base workspace:
%   sys, idxMV, idxMD, idxUD, mv_bounds, mu_u, sigma_u, mu_y, sigma_y,
%   x0, u0, p, options, Ts, varNames, outNames
%
% `options` (Mass/RelTol/AbsTol/MStateDependence) and `mv_bounds` are
% built the same way as in your Main_Exec.m, before the "Design MPC"
% section. run_nonlinear_test is assumed to have the exact signature
% you already use:
%   res = run_nonlinear_test(mpcobj, sys, p, options, idxMV, idxMD, ...
%           mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm, x0, u0, ...
%           r_traj, d_traj, Ts)
% returning res.y, res.u (and res.x if you use it elsewhere).
%
% This script does NOT modify design_mpc.m — it uses a local variant
% builder (build_mpc_variant) that mirrors it but exposes the knobs
% you want to sweep. Adjust the grid below to taste; a full-factorial
% sweep gets expensive fast (#configs x nonlinear-sim cost), so start
% small (8-16 configs) before widening it.
% ------------------------------------------------------------------

requiredVars = {'sys','idxMV','idxMD','idxUD','mv_bounds','mu_u','sigma_u', ...
    'mu_y','sigma_y','x0','u0','p','options','Ts'};
for kk = 1:numel(requiredVars)
    if ~evalin('base', sprintf('exist(''%s'',''var'')', requiredVars{kk}))
        error('sweep_mpc_settings:missingVar', ...
            'Required variable "%s" not found in base workspace. Run Main_Exec.m through the MPC design section first.', ...
            requiredVars{kk});
    end
    eval(sprintf('%s = evalin(''base'', ''%s'');', requiredVars{kk}, requiredVars{kk}));
end

%% ================= USER-CONFIGURABLE SWEEP GRID =================
NpList        = [30, 50, 75, 100, 150];                       % prediction horizons
NcList        = { [2 3 5 10 30 50], ...                        % your original blocking
                   [5 10 20 40], ...                            % coarser / fewer moves
                   [1 2 3 5 10 20 30 50], ...                   % finer near-term resolution
                   [10 30 60] };                                % very coarse, cheap to solve
rateScaleList = [0.25, 0.5, 1, 2, 4];                          % multiplies baseRateWeights
outWeightList = { [1 1], [2 1], [1 2], [5 1], [1 5], [3 1] };  % [T_weight, aH2O_weight]
rateLimitFracList = [0.005, 0.010, 0.020, 0.05];                % +/- fraction of MV range per sample
mvWeightScaleList = [0, 0.01, 0.1];                             % steady-state MV weight (0 = your original)

baseRateWeights = [10, 40, 10, 20, 1];  % your original ManipulatedVariablesRate

% Full factorial over the six dimensions above can easily reach
% thousands of configs (each needing its own nonlinear sim), so cap it:
% if the grid exceeds maxConfigs, randomly subsample instead of running
% everything. Set maxConfigs = Inf to force the full factorial.
maxConfigs = 40;
rngSeed = 1;   % fixed seed so the subsample is reproducible run-to-run

% Test scenario used to score every configuration: a moderate combined
% setpoint step (similar in spirit to S3), run on the NONLINEAR plant.
Tsim_test    = 500;
kStep_test   = 50;
deltaT_phys  = 1.0;      % K
deltaAH_phys = 0.005;    % water activity [-]
settleTolFrac = 0.02;    % settling band = 2% of step size

scoreWeights = struct( ...
    'IAE', 1.0, ...
    'Settle', 1.0, ...
    'Overshoot', 0.5, ...
    'Effort', 0.5, ...
    'SatPenalty', 5.0, ...   % added to normalized score if any MV saturates
    'FailPenalty', 1e6);     % score assigned to configs that error out

%% ================= BUILD PARAMETER GRID =================
paramSets = struct('Np',{},'Nc',{},'RateWeights',{},'OutputWeights',{}, ...
    'RateLimitFrac',{},'MVWeights',{},'label',{});
idxP = 0;
nu = numel(idxMV);  % number of MVs, pulled from your actual idxMV
for iNp = 1:numel(NpList)
    for iNc = 1:numel(NcList)
        for iRate = 1:numel(rateScaleList)
            for iOut = 1:numel(outWeightList)
                for iRL = 1:numel(rateLimitFracList)
                    for iMVw = 1:numel(mvWeightScaleList)
                        idxP = idxP + 1;
                        paramSets(idxP).Np = NpList(iNp);
                        paramSets(idxP).Nc = NcList{iNc};
                        paramSets(idxP).RateWeights = baseRateWeights * rateScaleList(iRate);
                        paramSets(idxP).OutputWeights = outWeightList{iOut};
                        paramSets(idxP).RateLimitFrac = rateLimitFracList(iRL);
                        paramSets(idxP).MVWeights = mvWeightScaleList(iMVw) * ones(1, nu);
                        paramSets(idxP).label = sprintf('Np%d_Nc%d_Rate%.2fx_Out%d-%d_RL%.3f_MVw%.2f', ...
                            NpList(iNp), numel(NcList{iNc}), rateScaleList(iRate), ...
                            outWeightList{iOut}(1), outWeightList{iOut}(2), ...
                            rateLimitFracList(iRL), mvWeightScaleList(iMVw));
                    end
                end
            end
        end
    end
end
fprintf('Full factorial grid size: %d configurations\n', numel(paramSets));

if numel(paramSets) > maxConfigs
    rng(rngSeed);
    keepIdx = sort(randperm(numel(paramSets), maxConfigs));
    paramSets = paramSets(keepIdx);
    fprintf(['Grid exceeds maxConfigs (%d) -> randomly subsampled to %d ' ...
        'configurations (seed=%d). Increase maxConfigs (or set it to Inf) ' ...
        'to run the full grid.\n'], maxConfigs, numel(paramSets), rngSeed);
end
fprintf('Total configurations to test: %d\n', numel(paramSets));

%% ================= REFERENCE TRAJECTORY FOR THE TEST SCENARIO =================
delta_norm_test = [deltaT_phys, deltaAH_phys] ./ sigma_y;
r_traj_test = zeros(Tsim_test, 2);
r_traj_test(kStep_test:end,:) = repmat(delta_norm_test, Tsim_test - kStep_test + 1, 1);
r_phys_test = r_traj_test .* sigma_y + mu_y;
t_vec_test  = (0:Tsim_test-1).' * Ts;

%% ================= SWEEP LOOP =================
results = table();
for k = 1:numel(paramSets)
    tune = paramSets(k);
    tStart = tic;
    try
        [mpc_k, mv_bounds_norm_k] = build_mpc_variant( ...
            sys, idxMV, idxMD, idxUD, mv_bounds, mu_u, sigma_u, Ts, tune);

        res = run_nonlinear_test(mpc_k, sys, p, options, idxMV, idxMD, ...
            mu_u, sigma_u, mu_y, sigma_y, mv_bounds_norm_k, x0, u0, ...
            r_traj_test, [], Ts);

        m = compute_metrics(t_vec_test, res.y, r_phys_test, res.u, mv_bounds, settleTolFrac);
        row = pack_result_row(tune, m, false, '');
    catch ME
        fprintf('Config %d/%d (%s) FAILED: %s\n', k, numel(paramSets), tune.label, ME.message);
        row = pack_result_row(tune, [], true, ME.message);
    end
    results = [results; row]; %#ok<AGROW>
    fprintf('[%d/%d] %-32s done in %.1fs\n', k, numel(paramSets), tune.label, toc(tStart));
end

%% ================= SCORE AND RANK =================
results = compute_scores(results, scoreWeights);
sortedResults = sortrows(results, 'Score');

fprintf('\n===== Top configurations (lower score = better) =====\n');
disp(sortedResults(1:min(10, height(sortedResults)), ...
    {'Label','Np','NcLen','RateScale','OutW_T','OutW_aH2O', ...
     'RateLimitFrac','MVWeight', ...
     'IAE_T','IAE_aH2O','Overshoot_T','Overshoot_aH2O', ...
     'SettleT_T','SettleT_aH2O','MVEffort','AnyMVSaturated','Score'}));

%% ================= PLOTS =================
nTop = min(10, height(sortedResults));
figure('Name','MPC Sweep: Top Configurations by Score');
barh(sortedResults.Score(nTop:-1:1));
set(gca, 'YTick', 1:nTop, 'YTickLabel', flipud(sortedResults.Label(1:nTop)), ...
    'TickLabelInterpreter','none');
xlabel('Composite score (lower is better)');
title(sprintf('Top %d of %d configurations', nTop, height(results)));
grid on;

valid = ~sortedResults.Failed;
figure('Name','MPC Sweep: Error vs Actuator Effort');
scatter(sortedResults.IAE_T(valid) + sortedResults.IAE_aH2O(valid), ...
    sortedResults.MVEffort(valid), 60, sortedResults.Np(valid), 'filled');
hold on; grid on;
xlabel('Combined IAE (T + a_{H2O})');
ylabel('Total MV effort (\Sigma|\Delta u|)');
title('Tracking error vs. actuator effort (color = N_p)');
cb = colorbar; cb.Label.String = 'Prediction horizon N_p';

% Best config to carry forward
bestTune = paramSets(strcmp({paramSets.label}, sortedResults.Label{1}));
fprintf('\nBest configuration: %s\n', sortedResults.Label{1});
disp(bestTune);

%% ================= LOCAL FUNCTIONS =================
function [mpcobj, mv_bounds_norm] = build_mpc_variant( ...
    sys, idxMV, idxMD, idxUD, mv_bounds, mu_u, sigma_u, Ts, tune)
% Mirrors design_mpc.m but exposes Np, Nc, RateWeights, OutputWeights
% as sweepable inputs. Keeps the same hard-constraint / scaling logic.

nu = numel(idxMV);

mv_bounds_norm = (mv_bounds - mu_u(idxMV).') ./ sigma_u(idxMV).';
range_u_norm = mv_bounds_norm(:,2) - mv_bounds_norm(:,1);

plant = setmpcsignals(sys, 'MV', idxMV, 'MD', idxMD, 'UD', idxUD);

mpcobj = mpc(plant, Ts, tune.Np, tune.Nc);

mpcobj.OV(1).ScaleFactor = 1;
mpcobj.OV(2).ScaleFactor = 1;

for i = 1:nu
    mpcobj.MV(i).ScaleFactor = range_u_norm(i);

    mpcobj.MV(i).Min = mv_bounds_norm(i,1);
    mpcobj.MV(i).Max = mv_bounds_norm(i,2);

    mpcobj.MV(i).RateMin = -tune.RateLimitFrac * range_u_norm(i);
    mpcobj.MV(i).RateMax =  tune.RateLimitFrac * range_u_norm(i);

    mpcobj.MV(i).MinECR = 0;
    mpcobj.MV(i).MaxECR = 0;

    mpcobj.MV(i).RateMinECR = 1;
    mpcobj.MV(i).RateMaxECR = 1;
end

mpcobj.Weights.OutputVariables = tune.OutputWeights;
mpcobj.Weights.ManipulatedVariables = tune.MVWeights;
mpcobj.Weights.ManipulatedVariablesRate = tune.RateWeights;

setEstimator(mpcobj, 'default');
end

function m = compute_metrics(t, y, r_phys, u, mv_bounds, settleTolFrac)
% All error/overshoot/settling metrics are computed per output channel
% (columns of y / r_phys). u and mv_bounds are per-MV (columns of u,
% rows of mv_bounds).

ny = size(y, 2);
dt = mean(diff(t));

m.IAE = zeros(1, ny);
m.overshoot = zeros(1, ny);
m.settlingTime = zeros(1, ny);

for kk = 1:ny
    e = r_phys(:,kk) - y(:,kk);
    m.IAE(kk) = sum(abs(e)) * dt;

    stepSize = r_phys(end,kk) - r_phys(1,kk);
    if abs(stepSize) > 1e-9
        if stepSize > 0
            m.overshoot(kk) = max(0, (max(y(:,kk)) - r_phys(end,kk)) / stepSize * 100);
        else
            m.overshoot(kk) = max(0, (r_phys(end,kk) - min(y(:,kk))) / stepSize * 100);
        end
    else
        m.overshoot(kk) = 0;
    end

    band = max(settleTolFrac * abs(stepSize), 1e-6);
    withinBand = abs(e) <= band;
    lastViolation = find(~withinBand, 1, 'last');
    if isempty(lastViolation)
        m.settlingTime(kk) = t(1);
    elseif lastViolation == numel(t)
        % Never settled within the simulated window: penalize rather
        % than propagate Inf into the scoring.
        m.settlingTime(kk) = 1.5 * t(end);
    else
        m.settlingTime(kk) = t(lastViolation + 1);
    end
end

m.MVeffort = sum(abs(diff(u, 1, 1)), 1);   % total variation per MV
nMV = size(u, 2);
m.MVsat = false(1, nMV);
for j = 1:nMV
    tol = 1e-6 * max(mv_bounds(j,2) - mv_bounds(j,1), eps);
    m.MVsat(j) = any(u(:,j) <= mv_bounds(j,1) + tol) || any(u(:,j) >= mv_bounds(j,2) - tol);
end
end

function row = pack_result_row(tune, m, failed, errMsg)
row = table();
row.Label = string(tune.label);
row.Np = tune.Np;
row.NcLen = numel(tune.Nc);
row.RateScale = tune.RateWeights(1) / 10;   % relative to base MV1 weight of 10
row.OutW_T = tune.OutputWeights(1);
row.OutW_aH2O = tune.OutputWeights(2);
row.RateLimitFrac = tune.RateLimitFrac;
row.MVWeight = tune.MVWeights(1);
row.Failed = failed;
row.ErrorMsg = string(errMsg);

if failed
    row.IAE_T = NaN; row.IAE_aH2O = NaN;
    row.Overshoot_T = NaN; row.Overshoot_aH2O = NaN;
    row.SettleT_T = NaN; row.SettleT_aH2O = NaN;
    row.MVEffort = NaN;
    row.AnyMVSaturated = 0;
else
    row.IAE_T = m.IAE(1);
    row.IAE_aH2O = m.IAE(2);
    row.Overshoot_T = m.overshoot(1);
    row.Overshoot_aH2O = m.overshoot(2);
    row.SettleT_T = m.settlingTime(1);
    row.SettleT_aH2O = m.settlingTime(2);
    row.MVEffort = sum(m.MVeffort);
    row.AnyMVSaturated = double(any(m.MVsat));
end
end

function results = compute_scores(results, w)
% Min-max normalize each raw metric across the VALID rows, combine into
% a single weighted score, and penalize failures / saturation heavily.

valid = ~results.Failed;
n = height(results);
results.Score = repmat(w.FailPenalty, n, 1);

if ~any(valid)
    warning('compute_scores:allFailed', 'Every configuration failed — nothing to score.');
    return;
end

combinedIAE    = results.IAE_T + results.IAE_aH2O;
combinedSettle = results.SettleT_T + results.SettleT_aH2O;
combinedOver   = results.Overshoot_T + results.Overshoot_aH2O;
effort         = results.MVEffort;

normalize = @(v) (v - min(v(valid))) ./ max(max(v(valid)) - min(v(valid)), eps);

score = w.IAE       * normalize(combinedIAE) ...
      + w.Settle    * normalize(combinedSettle) ...
      + w.Overshoot * normalize(combinedOver) ...
      + w.Effort    * normalize(effort) ...
      + w.SatPenalty * results.AnyMVSaturated;

results.Score(valid) = score(valid);
end