function [results, sysBest, best, diagData] = select_n4sid_model(zFull, varargin)
%SELECT_N4SID_MODEL Sweep N4SID orders and horizons using iddata input.
%
% Inputs
%   zFull : single-experiment iddata object containing normalized data
%
% Name-value options
%   'SplitFraction' : estimation-data fraction, default 0.70
%   'Hcand'          : candidate [r sy su] N4SID horizons
%   'OrderCand'      : candidate N4SID model orders
%   'OutputWeight'   : output weighting matrix for N4SID
%   'FitWeights'     : [w1 w2] weights on normalized Fit_y1/Fit_y2 when
%                      ranking, default [0.5 0.5]
%   'WriteCSV'       : save ranking table to CSV, default true
%   'CSVFile'        : CSV filename
%   'MakeFigures'    : generate compare/resid/IC figures, default true
%
% Outputs
%   results  : table containing all tested candidates
%   sysBest  : selected N4SID model
%   best     : one-row table containing selected candidate settings
%   diagData : struct containing zEst, zVal, compare options, sysAll, etc.

%% Input parsing
parser = inputParser;

defaultHcand = [
    5   10   10
    10   20   20
    10   30   30
    10   40   40
    20   30   30
    20   40   40
    20   60   60
    30   60   60
    30  100  100
    ];

addRequired(parser, 'zFull', @(z) isa(z, 'iddata'));
addParameter(parser, 'SplitFraction', 0.70, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(parser, 'Hcand', defaultHcand, ...
    @(x) isnumeric(x) && size(x,2) == 3 && all(x(:) > 0));
addParameter(parser, 'OrderCand', 6:2:24, ...
    @(x) isnumeric(x) && isvector(x) && all(x > 0));
addParameter(parser, 'OutputWeight', diag([0.3, 1]), ...
    @(x) isnumeric(x) && ismatrix(x));
addParameter(parser, 'FitWeights', [0.5, 0.5], ...
    @(x) isnumeric(x) && numel(x) == 2 && all(x >= 0));
addParameter(parser, 'WriteCSV', true, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'CSVFile', 'n4sid_horizon_order_results.csv', ...
    @(x) ischar(x) || isstring(x));
addParameter(parser, 'MakeFigures', true, @(x) islogical(x) && isscalar(x));

parse(parser, zFull, varargin{:});
cfg = parser.Results;

%% Ensure single-experiment data
if iscell(zFull.OutputData) || iscell(zFull.InputData)
    error(['select_n4sid_model currently supports one experiment only. ', ...
        'Pass one individual iddata experiment.']);
end

Ndata = size(zFull.OutputData, 1);
Ts = zFull.Ts;

if Ndata < 20
    error('The iddata object contains too few samples for N4SID validation.');
end

Nest = floor(cfg.SplitFraction * Ndata);
if Nest < 10 || Nest >= Ndata
    error('Invalid estimation/validation split.');
end

%% Split data chronologically
zEst = zFull(1:Nest);
zVal = zFull(Nest+1:end);

%% Candidate settings
Hcand = cfg.Hcand;
orderCand = cfg.OrderCand(:).';

% 'estimate' is not a valid compareOptions value -- use 'e'.
cmpValOpt = compareOptions('InitialCondition', 'e');
cmpZeroOpt = compareOptions('InitialCondition', 'z');
cmpEstOpt  = compareOptions('InitialCondition', 'e');

nTests = size(Hcand,1) * numel(orderCand);

idxCol = zeros(nTests,1);
rCol = zeros(nTests,1); syCol = zeros(nTests,1); suCol = zeros(nTests,1);
orderCol = zeros(nTests,1);
fitY1Col = nan(nTests,1); fitY2Col = nan(nTests,1);
stableCol = false(nTests,1);
controllableCol = false(nTests,1);
observableCol = false(nTests,1);
ctrbCondCol = nan(nTests,1);
obsvCondCol = nan(nTests,1);
negRealPoleCol = false(nTests,1);   % discrete pole near negative real axis
minTimeConstCol = nan(nTests,1);    % fastest continuous time constant [s]
maxTimeConstCol = nan(nTests,1);    % slowest continuous time constant [s]

sysAll = cell(nTests,1);   % cache every fitted model, avoid re-identifying later

%% N4SID sweep
q = 0;
for ih = 1:size(Hcand,1)
    for nx = orderCand
        q = q + 1;
        idxCol(q) = q;
        rCol(q) = Hcand(ih,1); syCol(q) = Hcand(ih,2); suCol(q) = Hcand(ih,3);
        orderCol(q) = nx;

        try
            optTest = n4sidOptions( ...
                'N4Weight', 'auto', ...
                'Focus', 'simulation', ...
                'N4Horizon', Hcand(ih,:), ...
                'EnforceStability', true, ...
                'InitialState', 'zero', ...
                'OutputWeight', cfg.OutputWeight);

            sysTest = n4sid(zEst, nx, optTest);
            sysAll{q} = sysTest;

            [~, fit] = compare(zVal, sysTest, cmpValOpt);
            fitY1Col(q) = fit(1);
            fitY2Col(q) = fit(2);

            A = sysTest.A; B = sysTest.B; C = sysTest.C;
            nxTest = size(A,1);

            stableCol(q) = all(abs(eig(A)) < 1);
            ctrbCondCol(q) = cond(ctrb(A,B));
            obsvCondCol(q) = cond(obsv(A,C));
            controllableCol(q) = rank(ctrb(A,B)) == nxTest;
            observableCol(q)   = rank(obsv(A,C)) == nxTest;

            % Discrete poles near the negative real axis indicate
            % near-Nyquist-frequency modes -- usually noise, not physics.
            z = eig(A);
            negRealPoleCol(q) = any(real(z) < 0 & abs(imag(z)) < 0.05*abs(real(z)+eps));

            % Continuous-time constants for physical-plausibility checks.
            % d2c warns verbosely when it hits real-negative discrete
            % poles (same condition flagged by negRealPoleCol above), so
            % suppress just that warning here to avoid spamming the
            % console once per candidate.
            ws = warning('off', 'all');
            try
                sysc = d2c(sysTest);
                Tc = -1 ./ real(eig(sysc.A));
                Tc = Tc(Tc > 0);
                if ~isempty(Tc)
                    minTimeConstCol(q) = min(Tc);
                    maxTimeConstCol(q) = max(Tc);
                end
            catch
                % d2c can fail outright on some candidates; leave NaN.
            end
            warning(ws);

        catch ME
            warning(['N4SID failed for H = [%d %d %d], order = %d.\n', ...
                'Reason: %s'], ...
                Hcand(ih,1), Hcand(ih,2), Hcand(ih,3), nx, ME.message);
        end
    end
end

%% Create results table
results = table( ...
    idxCol, rCol, syCol, suCol, orderCol, ...
    fitY1Col, fitY2Col, ...
    stableCol, controllableCol, observableCol, ...
    ctrbCondCol, obsvCondCol, negRealPoleCol, ...
    minTimeConstCol, maxTimeConstCol, ...
    'VariableNames', { ...
    'Idx', 'r', 'sy', 'su', 'Order', ...
    'Fit_y1', 'Fit_y2', ...
    'Stable', 'Controllable', 'Observable', ...
    'CtrbCond', 'ObsvCond', 'HasNegRealPole', ...
    'MinTimeConst_s', 'MaxTimeConst_s'});

results = results(~isnan(results.Fit_y2), :);
if isempty(results)
    error('All N4SID candidates failed. Check data length and horizon choices.');
end

% Combined score: normalize each fit column to [0,1] across candidates so
% that Fit_y1's large spread and Fit_y2's near-flat spread both contribute
% meaningfully, instead of the raw-percentage sort defaulting to whichever
% column happens to vary most across this dataset.
% NOTE: range() is computed manually below rather than calling the
% builtin, since a same-named function elsewhere on the path was
% shadowing MATLAB's range() and breaking this call.
f1 = results.Fit_y1; f2 = results.Fit_y2;
n1 = (f1 - min(f1)) / max(max(f1) - min(f1), eps);
n2 = (f2 - min(f2)) / max(max(f2) - min(f2), eps);
w = cfg.FitWeights / sum(cfg.FitWeights);
results.Score = w(1)*n1 + w(2)*n2;

% Hard filter: instability is disqualifying. Controllable/Observable are
% kept as informational columns (with condition numbers) rather than a
% second hard filter, since the binary rank test is unreliable at high
% order -- inspect CtrbCond/ObsvCond for borderline cases instead.
results = sortrows(results, {'Stable', 'Score'}, {'descend', 'descend'});

fprintf('\nTop N4SID candidates ranked by combined score (weights = [%.2f %.2f]):\n', w);
disp(results(1:min(15,height(results)), ...
    {'r','sy','su','Order','Fit_y1','Fit_y2','Score','Stable','HasNegRealPole','MinTimeConst_s'}));

%% Save results
if cfg.WriteCSV
    writetable(results, cfg.CSVFile);
    fprintf('N4SID ranking saved to: %s\n', cfg.CSVFile);
end

%% Selected candidate
best = results(1,:);

fprintf('\nSelected N4SID candidate:\n');
fprintf('N4Horizon = [%d %d %d]\n', best.r, best.sy, best.su);
fprintf('Order     = %d\n', best.Order);
fprintf('Fit y1    = %.2f %%\n', best.Fit_y1);
fprintf('Fit y2    = %.2f %%\n', best.Fit_y2);
fprintf('Score     = %.3f\n', best.Score);
fprintf('Stable    = %d | Controllable = %d | Observable = %d\n', ...
    best.Stable, best.Controllable, best.Observable);
if best.HasNegRealPole
    fprintf(2, 'WARNING: selected model has a discrete pole near the negative real axis.\n');
end

% Reuse the cached model instead of re-identifying (deterministic, so this
% would otherwise just repeat work already done in the loop above).
sysBest = sysAll{best.Idx};

%% Diagnostic figures
if cfg.MakeFigures

    figure('Name', 'Best N4SID Candidate: Validation');
    compare(zVal, sysBest, cmpValOpt);
    grid on;

    figure('Name', 'Best N4SID Candidate: Residual Analysis');
    resid(zVal, sysBest);
    grid on;

    figure('Name', 'Initial-Condition Diagnostic');

    subplot(2,1,1);
    compare(zFull, sysBest, cmpZeroOpt);
    grid on;
    title('Full data comparison: zero initial state');

    subplot(2,1,2);
    compare(zFull, sysBest, cmpEstOpt);
    grid on;
    title('Full data comparison: estimated initial state');

end

%% Return diagnostics
diagData = struct();
diagData.zFull = zFull;
diagData.zEst = zEst;
diagData.zVal = zVal;
diagData.Ts = Ts;
diagData.Ndata = Ndata;
diagData.Nest = Nest;
diagData.Nval = Ndata - Nest;
diagData.cmpValOpt = cmpValOpt;
diagData.cmpZeroOpt = cmpZeroOpt;
diagData.cmpEstOpt = cmpEstOpt;
diagData.config = cfg;
diagData.sysAll = sysAll;   % every candidate, in case you want to inspect another one

end