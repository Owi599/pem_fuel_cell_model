function plot_results(t, y_log, r_phys, u_log, mv_bounds, ...
    outNames, varNamesMV, figTitle, yHalfRange)
if nargin < 9 || isempty(yHalfRange)
    yHalfRange = [];
end
%% Output plots
figure('Name', figTitle);

for kk = 1:size(y_log,2)

    ax = subplot(size(y_log,2),1,kk);

    plot(t, y_log(:,kk), 'b', 'LineWidth', 1.5);
    hold on;
    grid on;

    plot(t, r_phys(:,kk), 'k--', 'LineWidth', 1.5);

    ax.FontSize = 22;
    ax.FontWeight = 'bold';
    ax.LineWidth = 1.0;



    % Apply only when a scenario provides engineering display limits.
    if ~isempty(yHalfRange)
        ylim(r_phys(1,kk) + [-yHalfRange(kk), yHalfRange(kk)]);
        % Limit y-axis tick labels to two decimal places.
        ytickformat('%.2f');
    end

    ylabel(outNames{kk}, ...
        'Interpreter', 'latex', ...
        'FontSize', 26, ...
        'FontWeight', 'bold');

    xlabel('Time [s]', ...
        'Interpreter', 'latex', ...
        'FontSize', 26, ...
        'FontWeight', 'bold');

    legend('Output', 'Reference', ...
        'Location', 'best', ...
        'Interpreter', 'latex', ...
        'FontSize', 22, ...
        'FontWeight', 'bold');
end
%% Manipulated-variable plots
if ~isempty(u_log)

    figure('Name', [figTitle ' — MVs']);

    for i = 1:size(u_log, 2)

        ax = subplot(size(u_log, 2), 1, i);

        plot(t, u_log(:,i), 'b', 'LineWidth', 1.4);
        hold on;
        grid on;

        yline(mv_bounds(i,1), 'r--', 'LineWidth', 1.2);
        yline(mv_bounds(i,2), 'r--', 'LineWidth', 1.2);

        % Tick-label appearance
        ax.FontSize = 13;
        ax.FontWeight = 'bold';
        ax.LineWidth = 1.0;

        ylabel(varNamesMV{i}, ...
            'Interpreter', 'none', ...
            'FontSize', 13, ...
            'FontWeight', 'bold');

        xlabel('Time [s]', ...
            'Interpreter', 'latex', ...
            'FontSize', 13, ...
            'FontWeight', 'bold');
    end
end

end