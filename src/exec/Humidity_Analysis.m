% Humidity Analysis around the Operating Point
% 0) Set Up
close all; clear; clc;
%%
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
targetI = [0.0162, 1, 2, 3, 4, 5,6,7,8];

%% 4) Compute steady-state simulation for each target current 
Tramp = 500; % [s] ramp duration for I_cell
Tss   = 5000; % [s] additional hold time to reach steady-state
Ttot  = Tramp + Tss;

nCols    = 60;   %  3 * p.N (anode, cathode, average humidity)
y        = zeros(numel(targetI),nCols);
x_ss_all = zeros(numel(targetI),numel(x0));

for i = 1:numel(targetI)
    u_target      = u0;
    u_target(iCellCol) = targetI(i);

    u_ramp = @(t) u0 + min(max(t,0),Tramp)/Tramp .* (u_target - u0);

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

    figure('Name', sprintf('Humidity, I_cell = %.3f A', targetI(i)));

    plot(zPos, aA, '-og', 'DisplayName', '$a^A_{H_2O}$','LineWidth',2); hold on;
    plot(zPos, aC, '-sb', 'DisplayName', '$a^C_{H_2O}$','LineWidth',2);
    plot(zPos, aAvg, '--k', 'DisplayName', '$\overline{a}_{H_2O}$','LineWidth',2.5);
   
    ax = gca;
    ax.FontSize = 20;
    ax.FontWeight= 'bold';
    %yline(1, ':k', 'saturation (a=1)', 'DisplayName', 'saturation limit');
    ylabel('Relative humidity [-]','FontSize',28,'FontWeight','bold');
   
    xlabel('Channel position z [mm]','FontSize',28,'FontWeight','bold');
    title(sprintf('Relative humidity (z) (%s) with I_{cell} = %.2f A', ...
        flowLabel, targetI(i)),'FontSize',28,'FontWeight','bold');
    legend('Location', 'best', 'Interpreter', 'latex','FontSize',24,'FontWeight','bold');
    grid on;
    fileName=sprintf('%s_I_cell_%.2f.svg',flowLabel,targetI(i));
        exportgraphics(ax,sprintf('%s',fileName),'Resolution',600)
end