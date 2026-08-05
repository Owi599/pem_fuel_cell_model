function plot_results(t, y_log, r_phys, u_log, mv_bounds, outNames, varNamesMV, figTitle)
figure('Name', figTitle);
for kk = 1:size(y_log,2)
    subplot(size(y_log,2),1,kk);
    plot(t, y_log(:,kk), 'b','LineWidth',1.2); hold on; grid on;
    plot(t, r_phys(:,kk), 'k--');
    ylabel(outNames{kk}, 'Interpreter','latex'); xlabel('Time [s]'); legend('Output','Reference');
end
if ~isempty(u_log)
    figure('Name',[figTitle ' — MVs']);
    for i = 1:size(u_log,2)
        subplot(size(u_log,2),1,i);
        plot(t, u_log(:,i)); hold on; grid on;
        yline(mv_bounds(i,1),'r--'); yline(mv_bounds(i,2),'r--');
        ylabel(varNamesMV{i});
    end
end
end