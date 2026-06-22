% Subspace ID data generation for PEMFC using ode15i/decic
close all; clear; clc;

%% 1) Parameters and options
p = mod_param_PEMFC();

options = odeset();
options.Mass = p.M();
options.RelTol = 1e-6;
options.AbsTol = 1e-6;
options.MStateDependence = 'none';

%% 2) Load model-input dataset
inputFileName   = 'dataset_100_251117_v06';
outputFileName_in = [inputFileName,'_inputs'];
rangeRows = [];

[u_traj, u_traj_info] = loadMatFile([outputFileName_in,'.mat'], ...
    rangeRows, p.inputs.variableNames);
uInterpolant_pp = griddedInterpolant(u_traj.time, u_traj.data, ...
    'pchip','nearest');

t = u_traj.time(:);
U = u_traj.data;
if size(U,1) ~= numel(t)
    U = U.';
end

Ts = mean(diff(t));
N  = size(U,1);
nu = size(U,2);

%% 3) Initial guess for state
% testbench initial input
in1 = [2.1, 560, 30, 70, 2.3, 164.73, 30, 70, 70, 0, 426, 2.7, 2.4];
initialInput2= p.testbench2struct(in1.');
x0 = steady_state_PEMFC(p, initialInput2,options);
u0 = U(1,:).';

%% 4) Input interpolant
u_fun = @(tt) interp1(t, U, tt, 'pchip', 'extrap').';
M = p.M();

%% 5) Define DAE residual: F(t,x,xdot)=M*xdot-f(x,u)=0
F = @(tt, x, xdot) pemfc_residual(tt, x, xdot, u_fun, p, M);

%% 6) Consistent initial conditions
x_guess    = x0(:);
xdot_guess = zeros(size(x_guess));

fixed_x    = false(size(x_guess));
fixed_xdot = false(size(x_guess));

[x0c, xdot0c] = decic(F, t(1), x_guess, fixed_x, xdot_guess, fixed_xdot, options);

%% 7) Simulate DAE with input trajectory
[t_sim, x_sim] = ode15i(F, t, x0c, xdot0c, options);

%% 8) Compute outputs
y_sim = sys_output_wrapper(x_sim, U, p);  

%% 9) Build iddata for N4SID
u_id = U;
y_id = y_sim;

u_id = u_id - mean(u_id(1:20,:), 1);
y_id = y_id - mean(y_id(1:20,:), 1);

z = iddata(y_id, u_id, Ts);

%% 10) Estimate linear model
nx = 4:20;   
sys = n4sid(z, 'best', 'N4Horizon', [10 10 10]);

compare(z, sys);