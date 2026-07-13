% Subspace ID data generation for PEMFC 
close all; clear; clc;
ensure_ode_pemfc_fresh();
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
y = zeros(N, 2);   

for k = 1:N
    y(k,:) = sys_output_wrapper(x(k,:).', U(k,:).', p).';
end
%% 8.1 Data Normalization
mu_u = mean(U,1); sigma_u = std(U,1);
mu_y = mean(y,1); sigma_y = std(y,1);

% guard against any zero-variance columns
sigma_u(sigma_u == 0) = 1;
sigma_y(sigma_y == 0) = 1;

normalize_u = @(U) (U - mu_u) ./ sigma_u;
normalize_y = @(Y) (Y - mu_y) ./ sigma_y;

U_Normalized = normalize_u(U);
y_Normalized = normalize_y(y);

fprintf('\nPer-channel training stats (mean, std):\n');
fprintf('  y1 (T_S(10))    : mean=%.3g, std=%.3g\n', mu_y(1), sigma_y(1));
fprintf('  y2 (a_H2O_avg)  : mean=%.3g, std=%.3g\n', mu_y(2), sigma_y(2));
%% 9) Estimated linear model

Ts = 1;

nx = 7:20;  
opt = n4sidOptions('Focus','simulation','N4Horizon',[30 30 30]);
sys = n4sid(U_Normalized,y_Normalized,nx,'Ts',Ts,opt);

A = sys.A;
B = sys.B; 
C = sys.C;
D = sys.D;

ev_sys = eig(A);
isStable = all(abs(ev_sys) < 1); % Stability condition for discrete-time systems 

% Display stability result
if isStable
    disp('The system is stable.');
else
    disp('The system is unstable.');
end
% Display controllability result

isControllable = rank(ctrb(A,B)) == min(size(ctrb(A,B)));
if isControllable
    disp('The system is controllable.');
else
    disp('The system is uncontrollable.');
end
%Display observability result

isObservable = rank(obsv(A,C)) == min(size(obsv(A,C)));
if isObservable
    disp('The system is observable.');
else
    disp('The system is unobservable.');
end

data_train = iddata(y,U,Ts);
data_train_Normalized = iddata(y_Normalized,U_Normalized,Ts);
%% 10) Validation Data set: different dataset, converted from testbench format
[u_traj_val, ~] = loadMatFile('dataset_50_OP_250622_v04_inputs.mat', [], ...
    p.testbench.variableNames);

u_traj_val_model = testbenchTraj2inputTraj(u_traj_val, p);

initialInput_val = p.testbench2struct(u_traj_val.data(1,:).');
x0_val = steady_state_PEMFC(p, initialInput_val, options);

uInterp_val = griddedInterpolant(u_traj_val_model.time, u_traj_val_model.data, ...
    'pchip', 'nearest');
[t_val, x_val] = ode15s(@(t,x) ode_PEMFC(t,x,uInterp_val(t).'), ...
    u_traj_val_model.time, x0_val, options);

U_val = u_traj_val_model.data;
if size(U_val,1) ~= numel(t_val); U_val = U_val.'; end

Nv = size(x_val,1);
y_val = zeros(Nv,2);
for k = 1:Nv
    y_val(k,:) = sys_output_wrapper(x_val(k,:).', U_val(k,:).', p).';
end

dt_val = mean(diff(t_val));
if abs(dt_val - Ts) > 1e-9
    warning('Held-out data sample spacing (%.3g s) differs from training (%.3g s).', dt_val, Ts);
end

U_val_d = U_val(1:Ts:end, :);
y_val_d = y_val(1:Ts:end, :);

% apply the SAME (training-derived) normalization
U_val_n = normalize_u(U_val_d);
y_val_n = normalize_y(y_val_d);

data_val = iddata(y_val_d, U_val_d, Ts);
data_val_n = iddata(y_val_n, U_val_n, Ts);


%% 11) Compare: training fit (baseline) vs. Validation test fit (the real test)/ Inf Compare Horizont 
figure('Name','Stage 0a (Inf Compare Horizont) - fit on TRAINING data (sanity baseline, expect high %)');
compare(data_train_Normalized, sys);

figure('Name','Stage 0b (Inf Compare Horizont) - fit on HELD-OUT data (the real generalization test)');
compare(data_val_n, sys,1);
[~, fit_val] = compare(data_val_n, sys);

fprintf('\nHeld-out fit (NRMSE %%) per output channel [T_S(10), a_H2O_avg]:\n');
disp(fit_val);

%% 12) Compare: training fit (baseline) vs. Validation test fit (the real test)/ 1-step Compare Horizont 
figure('Name','Stage 0a (1-step Compare Horizont) - fit on TRAINING data (sanity baseline, expect high %)');
compare(data_train_Normalized, sys);

figure('Name','Stage 0b (1-step Compare Horizont) - fit on HELD-OUT data (the real generalization test)');
compare(data_val_n, sys,1);
[~, fit_val_1_step] = compare(data_val_n, sys);

fprintf('\nHeld-out fit (NRMSE %%) per output channel [T_S(10), a_H2O_avg]:\n');
disp(fit_val_1_step);


%% LQI closed-loop simulation on one reduced model

nk = size(A,1);
ny = size(C,1);
nu = size(B,2);

Aaug = [A zeros(nk,ny);
    -C zeros(ny,ny)];
Baug = [B;
    -D];
Caug = [C zeros(ny,ny)];

isControllableaug = rank(ctrb(Aaug,Baug)) == min(size(ctrb(Aaug,Baug)));
if isControllableaug
    disp('The system is controllable.');
else
    disp('The system is uncontrollable.');
end

isObservableaug = rank(obsv(Aaug,Caug)) == min(size(obsv(Aaug,Caug)));
% Display observability result
if isObservableaug
    disp('The system is observable.');
else
    disp('The system is unobservable.');
end

Qx = 0.5*eye(nk);
Qi = 0.5*eye(ny);
Q  = blkdiag(Qx, Qi);
R  = 0.2*eye(nu);

Kaug = lqr(Aaug, Baug, Q, R);
Kx = Kaug(:,1:nk);
Ki = Kaug(:,nk+1:end);

%%
% Closed-loop simulation setup
Ts   = 100;
Tend = 100000;
t    = (0:Ts:Tend)';
Nsim = numel(t);

% Reference in deviation variables
r = zeros(Nsim, ny);
r(round(Nsim/3):end, :) = 1;   % step in both outputs, adjust as needed

x  = zeros(nk,1);
xi = zeros(ny,1);

X = zeros(nk, Nsim);
XI = zeros(ny, Nsim);
Y = zeros(ny, Nsim);
U = zeros(nu, Nsim);
E = zeros(ny, Nsim);

for k = 1:Nsim
    y = C*x ;
    e = r(k,:).' - y(:);

    u = -Kx*x - Ki*xi;

    X(:,k)  = x;
    XI(:,k) = xi;
    Y(:,k)  = y;
    U(:,k)  = u;
    E(:,k)  = e;

    xdot  = A*x + B*u;
    xidot = e;

    if k < Nsim
        x  = x  + Ts*xdot;
        xi = xi + Ts*xidot;
    end
end

figure;
subplot(3,1,1);
plot(t, Y.', 'LineWidth', 1.2);
grid on;
ylabel('y');
title('Closed-loop output tracking');

subplot(3,1,2);
plot(t, U.', 'LineWidth', 1.2);
grid on;
ylabel('u');

subplot(3,1,3);
plot(t, E.', 'LineWidth', 1.2);
grid on;
ylabel('e');
xlabel('Time [s]');