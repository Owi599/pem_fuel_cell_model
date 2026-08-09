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
%   'WriteCSV'       : save ranking table to CSV, default true
%   'CSVFile'        : CSV filename
%   'MakeFigures'    : generate compare/resid/IC figures, default true
%
% Outputs
%   results  : table containing all tested candidates
%   sysBest  : selected N4SID model
%   best     : one-row table containing selected candidate settings
%   diagData : struct containing zEst, zVal, compare options, etc.

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

addParameter(parser, 'WriteCSV', true, ...
    @(x) islogical(x) && isscalar(x));

addParameter(parser, 'CSVFile', 'n4sid_horizon_order_results.csv', ...
    @(x) ischar(x) || isstring(x));

addParameter(parser, 'MakeFigures', true, ...
    @(x) islogical(x) && isscalar(x));

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
% Estimation section begins at the known operating point.
zEst = zFull(1:Nest);

% Validation section is later in the same continuous trajectory.
zVal = zFull(Nest+1:end);

%% Candidate settings
Hcand = cfg.Hcand;
orderCand = cfg.OrderCand(:).';

cmpValOpt = compareOptions('InitialCondition', 'estimate');
cmpZeroOpt = compareOptions('InitialCondition', 'z');
cmpEstOpt = compareOptions('InitialCondition', 'estimate');

nTests = size(Hcand,1) * numel(orderCand);

rCol = zeros(nTests,1);
syCol = zeros(nTests,1);
suCol = zeros(nTests,1);
orderCol = zeros(nTests,1);

fitY1Col = nan(nTests,1);
fitY2Col = nan(nTests,1);
stableCol = false(nTests,1);
controllableCol = false(nTests,1);
observableCol = false(nTests,1);

q = 0;

%% N4SID sweep
for ih = 1:size(Hcand,1)

    for nx = orderCand
        q = q + 1;

        rCol(q) = Hcand(ih,1);
        syCol(q) = Hcand(ih,2);
        suCol(q) = Hcand(ih,3);
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

            [~, fit] = compare(zVal, sysTest, cmpValOpt);

            fitY1Col(q) = fit(1);
            fitY2Col(q) = fit(2);

            A = sysTest.A;
            B = sysTest.B;
            C = sysTest.C;
            nxTest = size(A,1);

            stableCol(q) = all(abs(eig(A)) < 1);
            controllableCol(q) = rank(ctrb(A,B)) == nxTest;
            observableCol(q) = rank(obsv(A,C)) == nxTest;

        catch ME
            warning(['N4SID failed for H = [%d %d %d], order = %d.\n', ...
                     'Reason: %s'], ...
                     Hcand(ih,1), Hcand(ih,2), Hcand(ih,3), nx, ME.message);
        end
    end
end

%% Create and rank results table
results = table( ...
    rCol, syCol, suCol, orderCol, ...
    fitY1Col, fitY2Col, ...
    stableCol, controllableCol, observableCol, ...
    'VariableNames', { ...
    'r', 'sy', 'su', 'Order', ...
    'Fit_y1', 'Fit_y2', ...
    'Stable', 'Controllable', 'Observable'});

% Remove failed identifications before ranking.
results = results(~isnan(results.Fit_y2), :);

if isempty(results)
    error('All N4SID candidates failed. Check data length and horizon choices.');
end

% Prefer stable, controllable, observable models; then maximize y2 and y1 fit.
results = sortrows(results, ...
    {'Stable', 'Controllable', 'Observable', 'Fit_y2', 'Fit_y1'}, ...
    {'descend', 'descend', 'descend', 'descend', 'descend'});

fprintf('\nTop N4SID candidates ranked by validation quality:\n');
disp(results(1:min(15,height(results)), :));

%% Save results
if cfg.WriteCSV
    writetable(results, cfg.CSVFile);
    fprintf('N4SID ranking saved to: %s\n', cfg.CSVFile);
end

%% Re-identify the selected model
best = results(1,:);

fprintf('\nSelected N4SID candidate:\n');
fprintf('N4Horizon = [%d %d %d]\n', best.r, best.sy, best.su);
fprintf('Order     = %d\n', best.Order);
fprintf('Fit y1    = %.2f %%\n', best.Fit_y1);
fprintf('Fit y2    = %.2f %%\n', best.Fit_y2);
fprintf('Stable    = %d\n', best.Stable);
fprintf('Controllable = %d\n', best.Controllable);
fprintf('Observable   = %d\n', best.Observable);

optBest = n4sidOptions( ...
    'N4Weight', 'auto', ...
    'Focus', 'simulation', ...
    'N4Horizon', [best.r, best.sy, best.su], ...
    'EnforceStability', true, ...
    'InitialState', 'zero', ...
    'OutputWeight', cfg.OutputWeight);

sysBest = n4sid(zEst, best.Order, optBest);

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

end