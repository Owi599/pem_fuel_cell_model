%N4SID selection base code
%% N4SID model-selection test: horizon and order sweep
Ndata = size(U_Normalized, 1);
Nest = floor(0.70 * Ndata);

% Estimation data starts at the operating point: zero initial state is valid.
zEst = iddata(y_Normalized(1:Nest,:), ...
              U_Normalized(1:Nest,:), Ts);

% Validation is the later part of the same continuous experiment.
% It generally does NOT start at the operating point.
zVal = iddata(y_Normalized(Nest+1:end,:), ...
              U_Normalized(Nest+1:end,:), Ts);

% [r sy su]: forward, past-output, and past-input horizons for N4SID.
Hcand = [
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

orderCand = 6:2:24;

% The validation segment starts away from the OP, so estimate its state.
cmpValOpt = compareOptions('InitialCondition', 'estimate');

nTests = size(Hcand,1) * numel(orderCand);

rCol      = zeros(nTests,1);
syCol     = zeros(nTests,1);
suCol     = zeros(nTests,1);
orderCol  = zeros(nTests,1);
fitY1Col  = nan(nTests,1);
fitY2Col  = nan(nTests,1);
stableCol = false(nTests,1);

q = 0;

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
                'OutputWeight', diag([0.3, 1]));

            sysTest = n4sid(zEst, nx, optTest);

            [~, fit] = compare(zVal, sysTest, cmpValOpt);

            fitY1Col(q) = fit(1);
            fitY2Col(q) = fit(2);
            stableCol(q) = all(abs(eig(sysTest.A)) < 1);

        catch ME
            warning('Failed: H=[%d %d %d], order=%d\n%s', ...
                Hcand(ih,1), Hcand(ih,2), Hcand(ih,3), nx, ME.message);
        end
    end
end

results = table(rCol, syCol, suCol, orderCol, fitY1Col, fitY2Col, stableCol, ...
    'VariableNames', {'r','sy','su','Order','Fit_y1','Fit_y2','Stable'});

% Highest y2 validation fit first; then use y1 to break ties.
results = sortrows(results, {'Fit_y2','Fit_y1'}, {'descend','descend'});

disp('Top N4SID candidates ranked by validation fit on y2:');
disp(results(1:min(15,height(results)), :));

writetable(results, 'n4sid_horizon_order_results.csv');

%% Rebuild and inspect the best identified candidate
best = results(1,:);

fprintf('\nBest validation candidate:\n');
fprintf('N4Horizon = [%d %d %d]\n', best.r, best.sy, best.su);
fprintf('Order     = %d\n', best.Order);
fprintf('Fit y1    = %.2f %%\n', best.Fit_y1);
fprintf('Fit y2    = %.2f %%\n', best.Fit_y2);

optBest = n4sidOptions( ...
    'N4Weight', 'auto', ...
    'Focus', 'simulation', ...
    'N4Horizon', [best.r best.sy best.su], ...
    'EnforceStability', true, ...
    'InitialState', 'zero', ...
    'OutputWeight', diag([0.3, 1]));

sysBest = n4sid(zEst, best.Order, optBest);

figure('Name', 'Best N4SID candidate on validation data');
compare(zVal, sysBest, cmpValOpt);
grid on;

figure('Name', 'Best N4SID candidate: residual analysis');
resid(zVal, sysBest);
grid on;
%% Initial-condition diagnostic on the full trajectory
zFull = iddata(y_Normalized, U_Normalized, Ts);

cmpZero = compareOptions('InitialCondition', 'z');
cmpEstimate = compareOptions('InitialCondition', 'estimate');

figure('Name', 'Initial-condition comparison on full data');

subplot(2,1,1);
compare(zFull, sysBest, cmpZero);
title('Comparison with zero initial state');

subplot(2,1,2);
compare(zFull, sysBest, cmpEstimate);
title('Comparison with estimated initial state');