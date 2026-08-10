% Humidity Analysis around the Operating Point
%% 0) Set Up
close all; clear; clc;
ensure_ode_pemfc_fresh();

%% 1) Parameters and options
p = mod_param_PEMFC();
options.Mass = p.M();
options.RelTol = 1e-6;
options.AbsTol = 1e-6;
options.MStateDependence = 'none';
flowLabel = 'co-flow';
if p.counter_flow
    flowLabel = 'counter flow';
end
disp(['Flow configuration: ', flowLabel]);

%% 2) Operating point from testbench initial input array
in1 = [2.1, 560, 30, 70, 2.3, 164.73, 30, 70, 70, 0, 426, 2.7, 2.4];
initialInput = p.testbench2struct(in1.');
x0 = steady_state_PEMFC(p, initialInput, options);
u0 = p.inputs2vec(testbench2model(initialInput, p));
u0 = u0(:);
iCellCol = find(strcmp(p.inputs.variableNames, 'I_cell'));

%% 3) Target current loads
targetI = [0.0162, 1, 2, 3, 4, 5, 6, 7, 8];

%% 4) Steady-state simulation sweep over I_cell
Tramp = 500;
Tss   = 5000;
Ttot  = Tramp + Tss;
nCols = 60;
y        = zeros(numel(targetI), nCols);
x_ss_all = zeros(numel(targetI), numel(x0));

for i = 1:numel(targetI)
    u_target = u0;
    u_target(iCellCol) = targetI(i);
    u_ramp = @(t) u0 + min(max(t,0),Tramp)/Tramp .* (u_target - u0);
    sol  = ode15s(@(t,x) ode_PEMFC(t,x,u_ramp(t)), [0 Ttot], x0, options);
    x_ss = deval(sol, Ttot);

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
zPos = (1:p.N) * p.delta_z * 1000;

%% 5) Output folder setup
folderpath = fullfile( ...
    'C:\Users\alsab\OneDrive\Desktop\UNI\Master\Master Thesis', ...
    'pem_fuel_cell_model', 'res', 'Figures', flowLabel);
if ~exist(folderpath, 'dir')
    mkdir(folderpath);
end

%% 6) Per-current profile plots
for i = 1:numel(targetI)
    aA   = y(i, 1:20);
    aC   = y(i, 21:40);
    aAvg = y(i, 41:60);
    figure('Name', sprintf('Humidity, I_cell = %.3f A', targetI(i)));
    plot(zPos, aA, '-og', 'DisplayName', '$a^A_{H_2O}$', 'LineWidth', 2); hold on;
    plot(zPos, aC, '-sb', 'DisplayName', '$a^C_{H_2O}$', 'LineWidth', 2);
    plot(zPos, aAvg, '--k', 'DisplayName', '$\overline{a}_{H_2O}$', 'LineWidth', 2.5);
    ax = gca; ax.FontSize = 20; ax.FontWeight = 'bold';
    ylabel('Relative humidity [-]', 'FontSize', 28, 'FontWeight', 'bold');
    xlabel('Channel position z [mm]', 'FontSize', 28, 'FontWeight', 'bold');
    title(sprintf('Relative humidity (z) (%s) with I_{cell} = %.2f A', ...
        flowLabel, targetI(i)), 'FontSize', 28, 'FontWeight', 'bold');
    legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
    grid on;
    fileName = sprintf('%s_I_cell_%.2f.svg', flowLabel, targetI(i));
    exportgraphics(gcf, fullfile(folderpath, fileName), 'ContentType', 'vector');
end

%% 7) Combined comparison plots across I_cell
cmap = turbo(numel(targetI));

% 7a) Anode
figure('Name', sprintf('Anode humidity comparison (%s)', flowLabel));
hold on;
for i = 1:numel(targetI)
    plot(zPos, y(i,1:20), '-o', 'Color', cmap(i,:), 'LineWidth', 2, ...
        'MarkerFaceColor', cmap(i,:), 'MarkerSize', 4);
end
ax = gca; ax.FontSize = 20; ax.FontWeight = 'bold';
ylabel('Anode relative humidity $a^A_{H_2O}$ [-]', 'FontSize', 26, 'FontWeight', 'bold', 'Interpreter', 'latex');
xlabel('Channel position z [mm]', 'FontSize', 26, 'FontWeight', 'bold');
title(sprintf('Anode humidity vs. channel position across $I_{cell}$ (%s)', flowLabel), ...
    'FontSize', 26, 'FontWeight', 'bold', 'Interpreter', 'latex');
grid on;
colormap(cmap); cb = colorbar;
cb.Label.String = 'I_{cell} [A]'; cb.Label.FontSize = 22; cb.Label.FontWeight = 'bold';
clim([1, numel(targetI)]); cb.Ticks = 1:numel(targetI);
cb.TickLabels = arrayfun(@(v) sprintf('%.2f', v), targetI, 'UniformOutput', false);
exportgraphics(gcf, fullfile(folderpath, sprintf('%s_anode_comparison_all_Icell.svg', flowLabel)), 'ContentType', 'vector');

% 7b) Cathode
figure('Name', sprintf('Cathode humidity comparison (%s)', flowLabel));
hold on;
for i = 1:numel(targetI)
    plot(zPos, y(i,21:40), '-s', 'Color', cmap(i,:), 'LineWidth', 2, ...
        'MarkerFaceColor', cmap(i,:), 'MarkerSize', 4);
end
ax = gca; ax.FontSize = 20; ax.FontWeight = 'bold';
ylabel('Cathode relative humidity $a^C_{H_2O}$ [-]', 'FontSize', 26, 'FontWeight', 'bold', 'Interpreter', 'latex');
xlabel('Channel position z [mm]', 'FontSize', 26, 'FontWeight', 'bold');
title(sprintf('Cathode humidity vs. channel position across $I_{cell}$ (%s)', flowLabel), ...
    'FontSize', 26, 'FontWeight', 'bold', 'Interpreter', 'latex');
grid on;
colormap(cmap); cb = colorbar;
cb.Label.String = 'I_{cell} [A]'; cb.Label.FontSize = 22; cb.Label.FontWeight = 'bold';
clim([1, numel(targetI)]); cb.Ticks = 1:numel(targetI);
cb.TickLabels = arrayfun(@(v) sprintf('%.2f', v), targetI, 'UniformOutput', false);
exportgraphics(gcf, fullfile(folderpath, sprintf('%s_cathode_comparison_all_Icell.svg', flowLabel)), 'ContentType', 'vector');

% 7c) Average
figure('Name', sprintf('Average humidity comparison (%s)', flowLabel));
hold on;
for i = 1:numel(targetI)
    plot(zPos, y(i,41:60), '--d', 'Color', cmap(i,:), 'LineWidth', 2, ...
        'MarkerFaceColor', cmap(i,:), 'MarkerSize', 4);
end
ax = gca; ax.FontSize = 20; ax.FontWeight = 'bold';
ylabel('Average relative humidity $\overline{a}_{H_2O}$ [-]', 'FontSize', 26, 'FontWeight', 'bold', 'Interpreter', 'latex');
xlabel('Channel position z [mm]', 'FontSize', 26, 'FontWeight', 'bold');
title(sprintf('Average humidity vs. channel position across $I_{cell}$ (%s)', flowLabel), ...
    'FontSize', 26, 'FontWeight', 'bold', 'Interpreter', 'latex');
grid on;
colormap(cmap); cb = colorbar;
cb.Label.String = 'I_{cell} [A]'; cb.Label.FontSize = 22; cb.Label.FontWeight = 'bold';
clim([1, numel(targetI)]); cb.Ticks = 1:numel(targetI);
cb.TickLabels = arrayfun(@(v) sprintf('%.2f', v), targetI, 'UniformOutput', false);
exportgraphics(gcf, fullfile(folderpath, sprintf('%s_average_comparison_all_Icell.svg', flowLabel)), 'ContentType', 'vector');

%% 8) Crossover analysis 
nI = numel(targetI);

% 8a) Anode–cathode crossover (per current, within-current comparison)
acTable = table();
for i = 1:nI
    [zc, vc] = find_crossovers(zPos, y(i,1:20), y(i,21:40));
    for k = 1:numel(zc)
        acTable = [acTable; table(targetI(i), zc(k), vc(k), ...
            'VariableNames', {'I_cell_A','z_mm','a_H2O'})];
    end
end
disp(acTable);
writetable(acTable, fullfile(folderpath, sprintf('%s_anode_cathode_crossovers.csv', flowLabel)));

% 8b) Anode–anode crossover (pairwise across all currents)
aaTable = table();
for i = 1:nI-1
    for j = i+1:nI
        [zc, vc] = find_crossovers(zPos, y(i,1:20), y(j,1:20));
        for k = 1:numel(zc)
            aaTable = [aaTable; table(targetI(i), targetI(j), zc(k), vc(k), ...
                'VariableNames', {'I1_A','I2_A','z_mm','a_H2O'})];
        end
    end
end
disp(aaTable);
writetable(aaTable, fullfile(folderpath, sprintf('%s_anode_anode_crossovers.csv', flowLabel)));

% 8c) Cathode–cathode crossover (pairwise across all currents)
ccTable = table();
for i = 1:nI-1
    for j = i+1:nI
        [zc, vc] = find_crossovers(zPos, y(i,21:40), y(j,21:40));
        for k = 1:numel(zc)
            ccTable = [ccTable; table(targetI(i), targetI(j), zc(k), vc(k), ...
                'VariableNames', {'I1_A','I2_A','z_mm','a_H2O'})];
        end
    end
end
disp(ccTable);
writetable(ccTable, fullfile(folderpath, sprintf('%s_cathode_cathode_crossovers.csv', flowLabel)));

%% 9) Peak-tracking analysis (anode)
peakVal = zeros(nI,1);
peakZ   = zeros(nI,1);
for i = 1:nI
    [peakVal(i), idxPeak] = max(y(i,1:20));
    peakZ(i) = zPos(idxPeak);
end
peakTable = table(targetI(:), peakVal, peakZ, ...
    'VariableNames', {'I_cell_A','PeakHumidity','PeakZ_mm'});
disp(peakTable);
writetable(peakTable, fullfile(folderpath, sprintf('%s_anode_peak_summary.csv', flowLabel)));

figure('Name', sprintf('Anode peak humidity vs I_cell (%s)', flowLabel));
plot(targetI, peakVal, '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.2 0.7 0.2], 'Color', [0.2 0.7 0.2]);
ax = gca; ax.FontSize = 20; ax.FontWeight = 'bold';
xlabel('$I_{cell}$ [A]', 'FontSize', 26, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('Peak anode humidity $\max_z\, a^A_{H_2O}$ [-]', 'FontSize', 26, 'FontWeight', 'bold', 'Interpreter', 'latex');
title(sprintf('Anode peak humidity vs. $I_{cell}$ (%s)', flowLabel), ...
    'FontSize', 26, 'FontWeight', 'bold', 'Interpreter', 'latex');
grid on;
exportgraphics(gcf, fullfile(folderpath, sprintf('%s_anode_peakvalue_vs_Icell.svg', flowLabel)), 'ContentType', 'vector');

figure('Name', sprintf('Anode peak location vs I_cell (%s)', flowLabel));
plot(targetI, peakZ, '-s', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.1 0.4 0.8], 'Color', [0.1 0.4 0.8]);
ax = gca; ax.FontSize = 20; ax.FontWeight = 'bold';
xlabel('$I_{cell}$ [A]', 'FontSize', 26, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('Peak location $z$ [mm]', 'FontSize', 26, 'FontWeight', 'bold');
title(sprintf('Anode peak humidity location vs. $I_{cell}$ (%s)', flowLabel), ...
    'FontSize', 26, 'FontWeight', 'bold', 'Interpreter', 'latex');
grid on;
exportgraphics(gcf, fullfile(folderpath, sprintf('%s_anode_peaklocation_vs_Icell.svg', flowLabel)), 'ContentType', 'vector');

figure('Name', sprintf('Anode peak summary (%s)', flowLabel));
yyaxis left
plot(targetI, peakVal, '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'Color', [0.2 0.7 0.2], 'MarkerFaceColor', [0.2 0.7 0.2]);
ylabel('Peak humidity [-]', 'FontSize', 24, 'FontWeight', 'bold');
ax = gca; ax.YColor = [0.2 0.7 0.2];
yyaxis right
plot(targetI, peakZ, '-s', 'LineWidth', 2, 'MarkerSize', 8, ...
    'Color', [0.1 0.4 0.8], 'MarkerFaceColor', [0.1 0.4 0.8]);
ylabel('Peak location z [mm]', 'FontSize', 24, 'FontWeight', 'bold');
ax = gca; ax.YColor = [0.1 0.4 0.8];
xlabel('$I_{cell}$ [A]', 'FontSize', 26, 'FontWeight', 'bold', 'Interpreter', 'latex');
title(sprintf('Anode peak humidity: value and location vs. $I_{cell}$ (%s)', flowLabel), ...
    'FontSize', 24, 'FontWeight', 'bold', 'Interpreter', 'latex');
ax = gca; ax.FontSize = 18; ax.FontWeight = 'bold';
grid on;
exportgraphics(gcf, fullfile(folderpath, sprintf('%s_anode_peak_combined_vs_Icell.svg', flowLabel)), 'ContentType', 'vector');

%% Helper Function
function [zCross, valCross] = find_crossovers(zPos, curve1, curve2)
    d = curve1(:) - curve2(:);
    signChange = find(diff(sign(d)) ~= 0);
    zCross = zeros(numel(signChange),1);
    valCross = zeros(numel(signChange),1);
    for k = 1:numel(signChange)
        i = signChange(k);
        z1 = zPos(i);   z2 = zPos(i+1);
        d1 = d(i);      d2 = d(i+1);
        frac = -d1 / (d2 - d1);
        zCross(k) = z1 + frac*(z2 - z1);
        v1 = curve1(i); v2 = curve1(i+1);
        valCross(k) = v1 + frac*(v2 - v1);
    end
end