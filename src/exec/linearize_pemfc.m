% Parameter set up
close all; clear; clc;

%% 1) Parameter and option set up
p = mod_param_PEMFC();

options = odeset();
options.Mass = p.M();
options.RelTol = 1e-6;
options.AbsTol = 1e-6;
options.MStateDependence = 'none';

%% 2) Load operating-point data

% inputFileName = 'dataset_15_OP_250622_v03';

inputFileName = 'dataset_100_251117_v06';
outputFileName_in = [inputFileName,'_inputs'];
rangeRows = [];

[u_traj_pp, u_traj_pp_info] = loadMatFile([outputFileName_in,'.mat'], ...
    rangeRows, p.inputs.variableNames);

% [u_traj, u_traj_info] = loadMatFile([outputFileName_in,'.mat'], ...
%     rangeRows, p.testbench.variableNames);

% u_traj_pp = testbenchTraj2inputTraj(u_traj, p);
uInterpolant_pp = griddedInterpolant(u_traj_pp.time, u_traj_pp.data, ...
    'pchip','nearest');

% testbench initial input
in1 = [2.1, 560, 30, 70, 2.3, 164.73, 30, 70, 70, 0, 426, 2.7, 2.4];

%Model Inputs
%initialInput = p.inputs2struct(u_traj_pp.data(1,:).');
%x0 = sys_states_PEMFC_initial(p, initialInput);
u0 = u_traj_pp.data(1,:).';

%Testbench Inputs
initialInput2= p.testbench2struct(in1.');
x0 = steady_state_PEMFC(p, initialInput2,options);
%u02 =  testbench2model(initialInput2,p);

%% Simulate nonlinear DAE
tspan = u_traj_pp.time;

u = @(t) uInterpolant_pp(t).';
[t, x] = ode15s(@(t,x) ode_PEMFC(t,x,u(t)), tspan, x0, options);

N = size(x,1);
y = zeros(N,2);
for k = 1:N
    y(k,:) = sys_output_wrapper(x(k,:).', u_traj_pp.data(k,:).', p).';
end
%% 3) Nonlinear model and dimensions
M  = full(p.M());
f  = @(x,u) ode_PEMFC(0,x,u);

nx = length(x0);
nu = length(u0);

A = zeros(nx,nx);
B = zeros(nx,nu);

f0 = f(x0,u0);
residual_inf = norm(f0,inf);

disp('--- Operating point check ---');
disp(['||f(x0,u0)||_inf = ', num2str(residual_inf)]);

%% 4) Numerical Jacobians with central differences and relative steps
sx = max(abs(x0),1);
su = max(abs(u0),1);

for i = 1:nx
    dx = zeros(nx,1);
    hx = 1e-6 * sx(i);
    dx(i) = hx;
    fp = f(x0 + dx, u0);
    fm = f(x0 - dx, u0);
    A(:,i) = (fp - fm) / (2*hx);
end

for j = 1:nu
    du = zeros(nu,1);
    hu = 1e-6 * su(j);
    du(j) = hu;
    fp = f(x0, u0 + du);
    fm = f(x0, u0 - du);
    B(:,j) = (fp - fm) / (2*hu);
end

E = M;
%% 5) Numerical Jacobians from the output function  for C and D matrices
g = @(x,u) sys_output_wrapper(x,u,p);
y0 = g(x0,u0);
size(x0)
size(x0.')
p.nstates
ny = length(y0);

C = zeros(ny,nx);
D = zeros(ny,nu);

sx = max(abs(x0),1);
su = max(abs(u0),1);

for i = 1:nx
    dx = zeros(nx,1);
    hx = 1e-6 * sx(i);
    dx(i) = hx;
    yp = g(x0 + dx, u0);
    ym = g(x0 - dx, u0);
    C(:,i) = (yp - ym)/(2*hx);
end

for j = 1:nu
    du = zeros(nu,1);
    hu = 1e-6 * su(j);
    du(j) = hu;
    yp = g(x0, u0 + du);
    ym = g(x0, u0 - du);
    D(:,j) = (yp - ym)/(2*hu);
end


%% 6) Sanity checks

% Perform sanity checks on the computed matrices
if any(isnan(A(:))) || any(isnan(B(:))) || any(isnan(C(:))) || any(isnan(D(:)))
    error('NaN values detected in the matrices.');
end

disp('--- Matrix sizes ---');
disp('Size of E:')
disp(size(E));
disp('Size of A:')
disp(size(A));
disp('Size of B:')
disp(size(B));
disp('Size of C:')
disp(size(C));
disp('Size of D:')
disp(size(D));

%% 7) Descriptor model and state space model.
 sysD= dss(A,B,C,D,E);
 sysE = dss2ss(sysD);
 figure;
 step(sysE)
 
 ev_sys = eig(sysE.A);
 
 isStable = all(real(ev_sys) < 0); 
 
 [A_obsv,B_obsv,C_obsv,T_obsv,k_obsv] = obsvf(sysE.A,sysE.B,sysE.C);

 rank_obsv = sum(k_obsv);
   
 [A_ctrb,B_ctrb,C_ctrb,T_ctrb,k_ctrb] = ctrbf(sysE.A,sysE.B,sysE.C);

 rank_ctrb = sum(k_ctrb);

     %% 12) Balanced truncation on explicit model only
    if isStable
        hsv = hsvd(sysE);

        figure;
        semilogy(hsv, 'o-');
        grid on;
        title('Hankel singular values');

        % simple automatic order suggestion: keep 99% HSV energy
        hsv_energy = cumsum(hsv) / sum(hsv);
        r_auto = find(hsv_energy >= 0.99, 1, 'first');
        r_auto = max(r_auto, 2);

        disp(['Suggested reduced order from HSV energy = ', num2str(r_auto)]);

        r_list = unique([10 15 r_auto]);
        r_list = r_list(r_list < order(sysE));

        redsys = cell(size(r_list));

        for k = 1:numel(r_list)
            redsys{k} = balred(sysE, r_list(k));
        end

        figure;
        step(sysE, redsys{:});
        grid on;
        leg = [{'full'}, arrayfun(@(r) ['r=',num2str(r)], r_list, 'UniformOutput', false)];
        legend(leg{:});
        title('Balanced truncation comparison');

        figure;
        bodemag(sysE, redsys{:});
        grid on;
        legend(leg{:});
        title('Frequency-response comparison');

        %% 13) Reduced-model diagnostics
        disp('--- Reduced model diagnostics ---');
        for k = 1:numel(r_list)
            rk = r_list(k);
            Ak = redsys{k}.A;
            Bk = redsys{k}.B;
            Ck = redsys{k}.C;
            Dk = redsys{k}.D;

            stable_k = all(real(eig(Ak)) < 0);
            [A_Co,B_Co,C_Co,T_Co,k_Co] = ctrbf(Ak,Bk,Ck);
            rankCo_k = sum(k_Co);
            [A_Ob,B_Ob,C_Ob,T_Ob,k_Ob] = ctrbf(Ak,Bk,Ck);
            rankOb_k = sum(k_Ob);

            disp(['r = ', num2str(rk), ...
                  ', stable = ', num2str(stable_k), ...
                  ', rank(ctrb) = ', num2str(rankCo_k), '/', num2str(size(Ak,1)), ...
                  ', rank(obsv) = ', num2str(rankOb_k), '/', num2str(size(Ak,1))]);
        end
    else
        warning('Explicit model is not asymptotically stable; skipping balanced truncation.');
    end
%%
sysr = redsys{2};
Ak = sysr.A;
Bk = sysr.B;
Ck = sysr.C;
Dk = sysr.D;
nk = size(Ak,1);
ny = size(Ck,1);
nu = size(Bk,2);
%%
% Linear vs nonlinear open-loop comparison

y_lin = lsim(sysr, u_traj_pp.data, t);
outNames = {'T_S(10) [K]','a_H2O_{avg} [-]'};
figure('Name','Linear (Jacobian) vs Nonlinear (DAE) open-loop response');
for kk = 1:2
    subplot(2,1,kk);
    plot(t, y(:,kk), 'g', 'LineWidth', 1.2); hold on;
    plot(t, y_lin(:,kk), 'r--', 'LineWidth', 1.2);
    ylabel(outNames{kk}); xlabel('Time [s]');
    legend('Nonlinear DAE','Linear (N4SID)');
end
figure('Name','compare(): NRMSE fit, linear model vs nonlinear DAE data');
compare(iddata(y, u_traj_pp.data, 1), sysr);


%% 12) LQI closed-loop simulation on one reduced model


% Augmented system for LQI
Aaug = [Ak zeros(nk,ny);
       -Ck zeros(ny,ny)];
Baug = [Bk;
       -Dk];
Caug = [Ck zeros(ny,ny)];

Qx = 0.5*eye(nk);
Qi = 0.5*eye(ny);
Q  = blkdiag(Qx, Qi);
R  = 0.2*eye(nu);

Kaug = lqr(Aaug, Baug, Q, R);
Kx = Kaug(:,1:nk);
Ki = Kaug(:,nk+1:end);

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
    y = Ck*x + Dk*zeros(nu,1);
    e = r(k,:).' - y(:);

    u = -Kx*x - Ki*xi;

    X(:,k)  = x;
    XI(:,k) = xi;
    Y(:,k)  = y;
    U(:,k)  = u;
    E(:,k)  = e;

    xdot  = Ak*x + Bk*u;
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

%% 13) MPC 

sysr = redsys{1};
Ak = sysr.A;
Bk = sysr.B;
Ck = sysr.C;
Dk = sysr.D;

G = ss(Ak,Bk,Ck,Dk);
G.D = zeros(size(G.D));
Ts = 100;
c = mpc(G,Ts);

c.Model.Nominal.Y = y0;
c.Model.Nominal.U = u0;
c.Model.Nominal.X = x0;

% Remove default output disturbance model.
setoutdist(c,"model",tf(0));   

% Set number of steps for simulation time of 5 seconds.
SimulationSteps = 1000;      

% Define Reference signal.
ref = y0+0.1*ones(SimulationSteps,1);

% Run closed-loop simulation.
[y, t] = sim(c,SimulationSteps,ref);


figure;
plot(t,ref,t,y);
xlabel('time');
title('Closed-Loop Response');
legend('ref','y');