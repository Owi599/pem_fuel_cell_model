% Humidity Analysis around the Operating Point
% 0) Set Up
close all; clear; clc;

ensure_ode_pemfc_fresh(); % making sure the mex execution file is up to 
                          % date to avoid overwriting custom parameters
                          % by the default

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

% Display the flow configuration
disp(['Flow configuration: ', flowLabel]);


%% 2) Computing Operating Point from testbench initail input array 

% testbench initial input
in1 = [2.1, 560, 30, 70, 2.3, 164.73, 30, 70, 70, 0, 426, 2.7, 2.4];

initialInput= p.testbench2struct(in1.');

x0 = steady_state_PEMFC(p, initialInput,options);
u0 = p.inputs2vec( testbench2model(initialInput, p) );

u0 = u0(:);

iCellCol = find(strcmp(p.inputs.variableNames, 'I_cell'));

%% 3) Target current loads
targetI = [0.01, 1, 2, 3, 4, 5, 6,7,8];

%% 4) Compute steady-state simulation for each target current 
Tramp = 2000; % [s] ramp duration for I_cell
Tss   = 4000; % [s] additional hold time to reach steady-state
Ttot  = Tramp + Tss;

nCols    = 80;   % 4 * p.N (anode, cathode, average, lambda_m)
y        = zeros(numel(targetI),nCols);
x_ss_all = zerso(numel(targetI),numel(x0));

for i = 1:numel(targetI)
    u_target      = u_op;
    u_target(iCellCol) = targetI(i);

    u_ramp = @(t) u_op + min(max(t,0),Tramp)/Tramp .* (u_target - u_op);

    sol  = ode15s(@(t,x) ode_PEMFC(t,x,u_ramp(t)), [0 Ttot], x0, options);
    x_ss = deval(sol, Ttot);

    % Convergence sanity check: state should barely move over the last
    % part of the hold if we have truly reached steady state.
    x_check   = deval(sol, Ttot - 500);
    relChange = max(abs(x_ss - x_check)) / max(1, max(abs(x_ss)));
    if relChange > 1e-3
        fprintf(['  Warning: I_cell = %.3f A steady state may not be fully ' ...
            'converged (rel. change over last 500 s = %.2e). Consider ' ...
            'increasing Tss.\n'], targetI(i), relChange);
    end

    x_ss_all(i,:) = x_ss.';
    y(i,:) = sys_output_wrapper_Analysis(x_ss, u_target, p).';
end

zPos = (1:p.N) * p.delta_z * 1000;   % physical channel position [mm]
%% 5) Visualization
for i = 1:numel(targetI)
    aA   = y(i, 1:20);
    aC   = y(i, 21:40);
    aAvg = y(i, 41:60);
    lam  = y(i, 61:80);

    figure('Name', sprintf('Humidity, I_cell = %.3f A', targetI(i)));

    yyaxis left;
    plot(zPos, aA, '-og', 'DisplayName', '$a^A_{H_2O}$'); hold on;
    plot(zPos, aC, '-sc', 'DisplayName', '$a^C_{H_2O}$');
    plot(zPos, aAvg, '--m', 'DisplayName', '$\overline{a}_{H_2O}$');
    yline(1, ':k', 'saturation (a=1)', 'DisplayName', 'saturation limit');
    ylabel('Water activity / relative humidity [-]');

    yyaxis right;
    plot(zPos, lam, '-.y', 'DisplayName', '$\lambda_m$', 'LineWidth', 1.3);
    ylabel('Membrane water content \lambda_m [-]');

    xlabel('Channel position z [mm]');
    title(sprintf('Along-channel humidity & membrane hydration (%s), I_{cell} = %.3f A -- fixed operating point, current-only sweep', ...
        flowLabel, targetI(i)));
    legend('Location', 'best', 'Interpreter', 'latex');
    grid on;
end