% Subspace ID data generation for PEMFC 
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
y = zeros(N, 2);   

for k = 1:N
    y(k,:) = sys_output_wrapper(x(k,:).', U(k,:).', p).';
end

%% 9) Time period

Ts = 100;


%% 10) Estimated linear model
nx = 7:20;   
sys = n4sid(U,y,nx,'Ts',Ts);

%%
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;
%%
ev_sys = eig(A);
isStable = all(abs(ev_sys) < 1); % Stability condition for discrete-time systems 

% Display stability result
if isStable
    disp('The system is stable.');
else
    disp('The system is unstable.');
end
%%

isControllable = rank(ctrb(A,B)) == min(size(ctrb(A,B)));
if isControllable
    disp('The system is controllable.');
else
    disp('The system is uncontrollable.');
end
%%

isObservable = rank(obsv(A,C)) == min(size(obsv(A,C)));
% Display observability result
if isObservable
    disp('The system is observable.');
else
    disp('The system is unobservable.');
end
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