function [y_meas, y_ctrl, yStruct, y_Analysis] = sys_output_PEMFC(x,u_traj_pp,p)
% --------------------------------------------------------------------------------------   
% model output function y = g(t,x,u)
% for assumed output vector (not test bench outputs!)
% y = [U_cell; T_a_out; T_c_out; T_S; p_a_out; p_c_out; ...
%           a_H2O_a_out; a_H2O_c_out]

% fnc inputs
    % x         - state vector
    % u_traj_pp - model input trajectories
    % p         - parameter structure
% fnc outputs
    % y         - output vector


% -------------------------------------------------------------------------
% output calculation from state vector and inputs
xStruct = p.state2struct(x);
if isstruct(u_traj_pp) && isfield(u_traj_pp,'data')
    u = p.inputs2struct(u_traj_pp.data.');
else 
    u = u_traj_pp; % Assign input trajectory directly if not a struct
end

% Cell voltage
yStruct.U_cell = xStruct.U_cell;

% output flows, temperature and pressure in the anode channel
if p.counter_flow

    idxA = 1;
    % counter-flow case
    yStruct.T_a_out = xStruct.T_a(idxA,:);
    yStruct.p_a_out = xStruct.p_a(p.N,:);
    yStruct.n_dot_H2_a_out = -2*p.K_a/p.delta_z*(u.p_a_out - xStruct.p_a(idxA,:)) .* xStruct.c_H2_a(idxA,:);
    yStruct.n_dot_H2O_a_out = -2*p.K_a/p.delta_z*(u.p_a_out - xStruct.p_a(idxA,:)) .* xStruct.c_H2O_a(idxA,:);
    yStruct.a_H2O_a_out = (p.R * xStruct.T_a(idxA,:) .* xStruct.c_H2O_a(idxA,:)) ./ ...
    p.p_sat(xStruct.T_a(idxA,:));
    yStruct.a_H2O_a_out_z = (p.R * xStruct.T_a .* xStruct.c_H2O_a) ./ ...
    p.p_sat(xStruct.T_a);
else
    
    idxA = p.N;
    % co-flow case
    yStruct.T_a_out = xStruct.T_a(idxA,:);
    yStruct.p_a_out = xStruct.p_a(1,:);
    yStruct.n_dot_H2_a_out = -2*p.K_a/p.delta_z*(u.p_a_out - xStruct.p_a(idxA,:)) .* xStruct.c_H2_a(idxA,:);
    yStruct.n_dot_H2O_a_out = -2*p.K_a/p.delta_z*(u.p_a_out - xStruct.p_a(idxA,:)) .* xStruct.c_H2O_a(idxA,:);
    yStruct.a_H2O_a_out = (p.R * xStruct.T_a(idxA,:) .* xStruct.c_H2O_a(idxA,:)) ./ ...
    p.p_sat(xStruct.T_a(idxA,:));
    yStruct.a_H2O_a_out_z = (p.R * xStruct.T_a .* xStruct.c_H2O_a) ./ ...
    p.p_sat(xStruct.T_a);
end

idxC = p.N;
% temperature and pressure in the cathode channel
yStruct.T_c_out = xStruct.T_c(idxC,:);
yStruct.p_c_out = xStruct.p_c(1,:);

% output flows in the cathode channel
yStruct.n_dot_O2_c_out = -2*p.K_c/p.delta_z*(u.p_c_out - xStruct.p_c(idxC,:)) .* xStruct.c_O2_c(idxC,:);
yStruct.n_dot_H2O_c_out = -2*p.K_c/p.delta_z*(u.p_c_out - xStruct.p_c(idxC,:)) .* xStruct.c_H2O_c(idxC,:);
yStruct.n_dot_N2_c_out = -2*p.K_c/p.delta_z*(u.p_c_out - xStruct.p_c(idxC,:)) .* xStruct.c_N2_c(idxC,:); % not used in the output vector

yStruct.a_H2O_c_out = (p.R * xStruct.T_c(idxC,:) .* xStruct.c_H2O_c(idxC,:)) ./ ...
p.p_sat(xStruct.T_c(idxC,:));
yStruct.a_H2O_c_out_z = (p.R * xStruct.T_c .* xStruct.c_H2O_c) ./ ...
    p.p_sat(xStruct.T_c);

% solid temperature/Humidity

yStruct.T_S = xStruct.T_s;
yStruct.Lambda = xStruct.lambda_m;

%Humidity
yStruct.c_H2O_c = xStruct.c_H2O_c;
yStruct.c_H2O_a = xStruct.c_H2O_a;
yStruct.a_H2O_avg = 0.5 * (yStruct.a_H2O_a_out + yStruct.a_H2O_c_out);

% --------------------------------------------------------------------------------------
% relative humidity at the two membrane interfaces (anode-side and cathode-side)
% mirrors the calculation in membrane.m, including the same smoothmax
% clipping applied there, so the reported values match what actually
% drives the membrane water flux rather than an unclipped surrogate
scaling = 300;
a_H2O_ca_raw = (xStruct.xi_H2O_ca .* xStruct.p_a) ./ p.p_sat(xStruct.T_s);
a_H2O_cc_raw = (xStruct.xi_H2O_cc .* xStruct.p_c) ./ p.p_sat(xStruct.T_s);

yStruct.a_H2O_ca_z = -smoothmax(-a_H2O_ca_raw, -0.97, scaling);
yStruct.a_H2O_cc_z = -smoothmax(-a_H2O_cc_raw, -0.97, scaling);

yStruct.a_H2O_avg_z = 0.5 * (yStruct.a_H2O_a_out_z + yStruct.a_H2O_c_out_z);

y_meas = [ ...
    yStruct.U_cell;
    yStruct.T_a_out;
    yStruct.T_c_out;
    yStruct.T_S;
    yStruct.p_a_out;
    yStruct.p_c_out;
    yStruct.a_H2O_a_out;
    yStruct.a_H2O_c_out ];

y_ctrl = [yStruct.T_S(10); 
    yStruct.a_H2O_avg 
    ];

y_Analysis = [yStruct.a_H2O_a_out_z;
    yStruct.a_H2O_c_out_z;
    yStruct.a_H2O_avg_z;
    yStruct.a_H2O_ca_z;
    yStruct.a_H2O_cc_z
    ];
end

