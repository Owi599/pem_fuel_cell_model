% Humidity Analysis around the Operating Point
% 0) Set Up
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
%% 5 ) Visualization Set up
% Output folder for the current flow case
    folderpath = fullfile( ...
        'C:\Users\alsab\OneDrive\Desktop\UNI\Master\Master Thesis', ...
        'pem_fuel_cell_model', 'res', 'Figures', flowLabel);

    % Create folder if it does not yet exist
    if ~exist(folderpath, 'dir')
        mkdir(folderpath);
    end
%% 6) Visualization
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
    

    % Keep the desired file naming scheme
    fileName = sprintf('%s_I_cell_%.2f.svg', flowLabel, targetI(i));
    filepath = fullfile(folderpath, fileName);

    % Export to the specified folder
    exportgraphics(gcf, filepath, 'ContentType', 'vector');
end
%% 7) Combined anode humidity comparison across I_cell (color-coded)
figure('Name', sprintf('Anode humidity comparison (%s)', flowLabel));
cmap = turbo(numel(targetI));   % perceptually-uniform colormap, one color per current
hold on;
for i = 1:numel(targetI)
    aA = y(i, 1:20);
    plot(zPos, aA, '-o', 'Color', cmap(i,:), 'LineWidth', 2, ...
        'MarkerFaceColor', cmap(i,:), 'MarkerSize', 4);
end
ax = gca;
ax.FontSize = 20;
ax.FontWeight = 'bold';
ylabel('Anode relative humidity $a^A_{H_2O}$ [-]', 'FontSize', 26, ...
    'FontWeight', 'bold', 'Interpreter', 'latex');
xlabel('Channel position z [mm]', 'FontSize', 26, 'FontWeight', 'bold');
title(sprintf('Anode humidity vs. channel position across $I_{cell}$ (%s)', flowLabel), ...
    'FontSize', 26, 'FontWeight', 'bold', 'Interpreter', 'latex');
grid on;

colormap(cmap);
cb = colorbar;
cb.Label.String = 'I_{cell} [A]';
cb.Label.FontSize = 22;
cb.Label.FontWeight = 'bold';
% Map colorbar ticks back to actual target currents rather than 1:N index
clim([1, numel(targetI)]);
cb.Ticks = 1:numel(targetI);
cb.TickLabels = arrayfun(@(v) sprintf('%.2f', v), targetI, 'UniformOutput', false);

fileNameAnode = sprintf('%s_anode_comparison_all_Icell.svg', flowLabel);
filepathAnode = fullfile(folderpath, fileNameAnode);

exportgraphics(gcf, filepathAnode, 'ContentType', 'vector');
%% 7b) Combined cathode humidity comparison across I_cell (color-coded)
figure('Name', sprintf('Cathode humidity comparison (%s)', flowLabel));
hold on;
for i = 1:numel(targetI)
    aC = y(i, 21:40);
    plot(zPos, aC, '-s', 'Color', cmap(i,:), 'LineWidth', 2, ...
        'MarkerFaceColor', cmap(i,:), 'MarkerSize', 4);
end
ax = gca;
ax.FontSize = 20;
ax.FontWeight = 'bold';
ylabel('Cathode relative humidity $a^C_{H_2O}$ [-]', 'FontSize', 26, ...
    'FontWeight', 'bold', 'Interpreter', 'latex');
xlabel('Channel position z [mm]', 'FontSize', 26, 'FontWeight', 'bold');
title(sprintf('Cathode humidity vs. channel position across $I_{cell}$ (%s)', flowLabel), ...
    'FontSize', 26, 'FontWeight', 'bold', 'Interpreter', 'latex');
grid on;

colormap(cmap);
cb = colorbar;
cb.Label.String = 'I_{cell} [A]';
cb.Label.FontSize = 22;
cb.Label.FontWeight = 'bold';
clim([1, numel(targetI)]);
cb.Ticks = 1:numel(targetI);
cb.TickLabels = arrayfun(@(v) sprintf('%.2f', v), targetI, 'UniformOutput', false);

fileNameCathode = sprintf('%s_cathode_comparison_all_Icell.svg', flowLabel);
filepathCathode = fullfile(folderpath, fileNameCathode);

exportgraphics(gcf, filepathCathode, 'ContentType', 'vector');
%% 7c) Combined average humidity comparison across I_cell (color-coded)
figure('Name', sprintf('Average humidity comparison (%s)', flowLabel));
hold on;
for i = 1:numel(targetI)
    aAvg = y(i, 41:60);
    plot(zPos, aAvg, '--d', 'Color', cmap(i,:), 'LineWidth', 2, ...
        'MarkerFaceColor', cmap(i,:), 'MarkerSize', 4);
end
ax = gca;
ax.FontSize = 20;
ax.FontWeight = 'bold';
ylabel('Average relative humidity $\overline{a}_{H_2O}$ [-]', 'FontSize', 26, ...
    'FontWeight', 'bold', 'Interpreter', 'latex');
xlabel('Channel position z [mm]', 'FontSize', 26, 'FontWeight', 'bold');
title(sprintf('Average humidity vs. channel position across $I_{cell}$ (%s)', flowLabel), ...
    'FontSize', 26, 'FontWeight', 'bold', 'Interpreter', 'latex');
grid on;

colormap(cmap);
cb = colorbar;
cb.Label.String = 'I_{cell} [A]';
cb.Label.FontSize = 22;
cb.Label.FontWeight = 'bold';
clim([1, numel(targetI)]);
cb.Ticks = 1:numel(targetI);
cb.TickLabels = arrayfun(@(v) sprintf('%.2f', v), targetI, 'UniformOutput', false);

fileNameAvg = sprintf('%s_average_comparison_all_Icell.svg', flowLabel);
filepathAvg = fullfile(folderpath, fileNameAvg);

exportgraphics(gcf, filepathAvg, 'ContentType', 'vector');