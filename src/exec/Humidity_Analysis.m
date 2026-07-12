% Humidity Analysis
close all; clear; clc;

%% 1) Parameters and options
p = mod_param_PEMFC();

options.Mass=p.M();
options.RelTol=1e-6;
options.AbsTol=1e-6;
options.MStateDependence='none';

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


%% 7) Simulate DAE with input trajectory
[t, x] = ode15s(@(t,x) ode_PEMFC(t,x,u(t)), tspan, x0, options);

%% 8) Compute outputs
N = size(x,1);
y = zeros(N, 60);   

for k = 1:N
    y(k,:) = sys_output_wrapper_Analysis(x(k,:).', U(k,:).', p).';
end
%% I_cell = 0.0162 A
figure;
plot(y(17282,1:20),'DisplayName','$a^A_{H2O}$')
hold on;
plot(y(17282,21:40),'DisplayName','$a^C_{H2O}$')
plot(y(17282,41:60),'DisplayName','$\overline{a}_{H20}$')
xlim([1,20]);
xlabel('Z');
ylabel('Relative Humidity');
title('Relative Humidity over different Current Loads');
legend('Location','best',Interpreter='latex');
grid on;

%% I_cell = 1.0003 A
figure;
plot(y(9201,1:20),'DisplayName','$a^A_{H2O}$')
hold on;
xlim([1,20]);
plot(y(9201,21:40),'DisplayName','$a^C_{H2O}$')
plot(y(9201,41:60),'DisplayName','$\overline{a}_{H20}$')
xlabel('Z');
ylabel('Relative Humidity');
title('Relative Humidity over different Current Loads');
legend('Location','best',Interpreter='latex');
grid on;
%% I_cell = 2.0939 A
figure;
plot(y(8875,1:20),'DisplayName','$a^A_{H2O}$')
hold on;
xlim([1,20]);
plot(y(8875,21:40),'DisplayName','$a^C_{H2O}$')
plot(y(8875,41:60),'DisplayName','$\overline{a}_{H20}$')
xlabel('Z');
ylabel('Relative Humidity');
title('Relative Humidity over different Current Loads');
legend('Location','best',Interpreter='latex');
grid on;
%% I_cell = 3.0006 A
figure;
plot(y(15691,1:20),'DisplayName','$a^A_{H2O}$')
hold on;
xlim([1,20]);
plot(y(15691,21:40),'DisplayName','$a^C_{H2O}$')
plot(y(15691,41:60),'DisplayName','$\overline{a}_{H20}$')
xlabel('Z');
ylabel('Relative Humidity');
title('Relative Humidity over different Current Loads');
legend('Location','best',Interpreter='latex');
grid on;
%% I_cell = 3.9928 A
figure;
plot(y(4885,1:20),'DisplayName','$a^A_{H2O}$')
hold on;
xlim([1,20]);
plot(y(4885,21:40),'DisplayName','$a^C_{H2O}$')
plot(y(4885,41:60),'DisplayName','$\overline{a}_{H20}$')
xlabel('Z');
ylabel('Relative Humidity');
title('Relative Humidity over different Current Loads');
legend('Location','best',Interpreter='latex');
grid on;
%% I_cell = 4.9946 A
figure;
plot(y(4541,1:20),'DisplayName','$a^A_{H2O}$')
hold on;
xlim([1,20]);
plot(y(4541,21:40),'DisplayName','$a^C_{H2O}$')
plot(y(4541,41:60),'DisplayName','$\overline{a}_{H20}$')
xlabel('Z');
ylabel('Relative Humidity');
title('Relative Humidity over different Current Loads');
legend('Location','best',Interpreter='latex');
grid on;
%% I_cell = 5.7628 A
figure;
plot(y(11021,1:20),'DisplayName','$a^A_{H2O}$')
hold on;
xlim([1,20]);
plot(y(11021,21:40),'DisplayName','$a^C_{H2O}$')
plot(y(11021,41:60),'DisplayName','$\overline{a}_{H20}$')
xlabel('Z');
ylabel('Relative Humidity');
title('Relative Humidity over different Current Loads');
legend('Location','best',Interpreter='latex');
grid on;

