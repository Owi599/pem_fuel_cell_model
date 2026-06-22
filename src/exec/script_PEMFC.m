%%%
%%% reference fuel cell model %%%%%%%%%%%%%%
%%%

%   - fuel cell modeling according to [Mangold2010]
%   - the model is discretized along the channel
%   - formation of liquid water is not considered
%   - nonlinear DAE system
%   - initial conditions from function
%   - model inputs can be defined as function of time and/or state
%   - JPattern and sparse mass matrix for efficient computations
%   - possible coupling of test bench data to model
% !!! relative humidity at membrane limited to <1 in mod_param_PEMFC.m
% !!! empiric membrane correlations corrected for low humidity

%% initialization
clear all; close all; clc;

% use matlab coder for faster computation by generating mex-function of ode PEMFC

force_codegen = 0; % 0 - use presaved mex-file; 1 - new mex-file needed

if force_codegen 
    ode_PEMFC_CoderScript();
end

% model parameters
p = mod_param_PEMFC();

%% options

% options
force_rec_jac = 0; % 0 - use presaved JPattern; 1 - new JPattern needed
dequant_IO = 0; % 1 - use de-quantization filter for input trajectories; 0 - don't use it

% options for ode15s-solver
options.Mass=p.M();
options.RelTol=1e-6;
options.AbsTol=1e-6;
options.MStateDependence='none';

JPattern = myJPattern();
rec_jac = force_rec_jac || p.nstates ~= size(JPattern,1);

if rec_jac == 0
    options.JPattern = JPattern;
end

%% load input trajectory data

inputFileName = 'dataset_100_251117_v06';
outputFileName_in = [inputFileName,'_inputs'];

% manually select a specific part of the experiment (the selected section is assumed to start in a steady state)
rangeRows = [];% row indizes for data, '[]' for entire experiment; 

% load testbench inputs
%[u_traj,u_traj_info] = loadMatFile([outputFileName_in,'.mat'],rangeRows,p.testbench.variableNames);
[u_traj,u_traj_info] = loadMatFile([outputFileName_in,'.mat'],rangeRows,p.inputs.variableNames);

%% pre-processing and interpolation for data
info('Pre-processing started.');

% de-quantize input trajectories
if dequant_IO
    for i1=1:length(u_traj_info.variableNames)
        [u_traj.data(:,i1),~] = deQuant(u_traj.data(:,i1));
    end
end

% convert testbench inputs to model inputs
u_traj_pp = testbenchTraj2inputTraj(u_traj, p);

% create interpolant for model inputs
uInterpolant_pp = griddedInterpolant(u_traj_pp.time, u_traj_pp.data,'pchip','nearest'); % spline, pchip, linear

info('Pre-processing complete.');

%% state vector initialization
% calculate steady state from initial inputs
info('Initialization started.');

% initial point from measurement file 
initialInput = p.testbench2struct(u_traj.data(1,:).');

% initial conditions for experiments
x0_initialization = steady_state_PEMFC(p,initialInput,options);

info('Initialization complete.');

%% JPattern
% generate JPattern if it is needed
if rec_jac == 1
    info('JPattern calculation started.', 1);
    
    % calculate model input vector from testbench input struct
    u_vec = p.inputs2vec(testbench2model(initialInput,p));
    
    % generate and save new Jpattern
    genJpattern(@(x) ode_PEMFC(0,x,u_vec),x0_initialization);
    options.JPattern=myJPattern();
    
    info('JPattern calculation complete.', 1);
end


%% experiment

info('Simulation started.');

% time span [s]
tspan = u_traj.time;

% initial value from initialization
x0 = x0_initialization;

% innput values from interpolant
u = @(t) uInterpolant_pp(t).';


% solution with ode15s-solver
[t, x] = ode15s(@(t,x) ode_PEMFC(t,x,u(t)), tspan, x0, options);

info('Simulation complete.');

% output calculation
y = sys_output_PEMFC(x,u_traj_pp,p);

%% plots
% voli - volume number
for voli = [1,round(p.N/2),p.N]
    plots_PEMFC_0D_OPcheck(x,t,u_traj_pp,voli,p);
end
