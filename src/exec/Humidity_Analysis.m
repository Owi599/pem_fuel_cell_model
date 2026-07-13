% Humidity Analysis
close all; clear; clc;
ensure_ode_pemfc_fresh();

%% 1) Parameters and options
p = mod_param_PEMFC();

options.Mass=p.M();
options.RelTol=1e-6;
options.AbsTol=1e-6;
options.MStateDependence='none';

flowLabel = 'co-flow';

if p.counter_flow

    flowLabel='counter flow';

end
disp(p.counter_flow)

%% 2) Load model-input dataset
inputFileName   = 'dataset_100_251117_v06';
outputFileName_in = [inputFileName,'_inputs'];
rangeRows = [];

[u_traj, u_traj_info] = loadMatFile([outputFileName_in,'.mat'], ...
    rangeRows, p.inputs.variableNames);
uInterpolant_pp = griddedInterpolant(u_traj.time, u_traj.data, ...
    'pchip','nearest');
%% 3) Compute Initial steady state
% testbench initial input
in1 = [2.1, 560, 30, 70, 2.3, 164.73, 30, 70, 70, 0, 426, 2.7, 2.4];
initialInput= p.testbench2struct(in1.');
x0 = steady_state_PEMFC(p, initialInput,options);
u0 = u_traj.data(1,:).';
tspan = u_traj.time;
U = u_traj.data;
if size(U,1) ~= numel(tspan)
    U = U.';
end
%% 4) Input values from interpolent
u = @(t) uInterpolant_pp(t).';


%% 5) Simulate DAE with input trajectory
[t, x] = ode15s(@(t,x) ode_PEMFC(t,x,u(t)), tspan, x0, options);
%% 6) Locate the I_cell column in the input structure
iCellCol = find(strcmp(p.inputs.variableNames, 'I_cell'));
%% 7) Target current loads -- indices found dynamically, one per target
targetI = [0, 1, 2, 3,3.54, 4, 5, 6];
idxList = zeros(size(targetI));
for i = 1:numel(targetI)
    [~, idxList(i)] = min(abs(U(:,iCellCol) - targetI(i)));
end
%% 8) Compute outputs
nCols = 80;   % 4 * p.N (anode, cathode, average, lambda_m)
y = zeros(numel(idxList), nCols);
for i = 1:numel(idxList)
    k = idxList(i);
    y(i,:) = sys_output_wrapper_Analysis(x(k,:).', U(k,:).', p).';
end

zPos = (1:p.N) * p.delta_z * 1000;   % physical channel position [mm]
%% 9) Visualization
for i = 1:numel(targetI)
    aA   = y(i, 1:20);
    aC   = y(i, 21:40);
    aAvg = y(i, 41:60);
    lam  = y(i, 61:80);

    I_actual = U(idxList(i), iCellCol);

    figure('Name', sprintf('Humidity, I_cell = %.3f A', I_actual));

    yyaxis left;
    plot(zPos, aA, '-og', 'DisplayName', '$a^A_{H_2O}$'); hold on;
    plot(zPos, aC, '-sm', 'DisplayName', '$a^C_{H_2O}$');
    plot(zPos, aAvg, '--c', 'DisplayName', '$\overline{a}_{H_2O}$');
    yline(1, ':k', 'saturation (a=1)', 'DisplayName', 'saturation limit');
    ylabel('Water activity / relative humidity [-]');

    yyaxis right;
    plot(zPos, lam, '-.', 'DisplayName', '$\lambda_m$', 'LineWidth', 1.3);
    ylabel('Membrane water content \lambda_m [-]');

    xlabel('Channel position z [mm]');
    title(sprintf('Along-channel humidity & membrane hydration (%s), I_{cell} = %.3f A', ...
        flowLabel, I_actual));
    legend('Location', 'best', 'Interpreter', 'latex');
    grid on;
end