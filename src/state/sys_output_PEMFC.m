function [y_meas, y_ctrl, yStruct] = sys_output_PEMFC(x,u_traj_pp,p)
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

else
    
    idxA = p.N;
    % co-flow case
    yStruct.T_a_out = xStruct.T_a(idxA,:);
    yStruct.p_a_out = xStruct.p_a(p.N,:);
    yStruct.n_dot_H2_a_out = -2*p.K_a/p.delta_z*(u.p_a_out - xStruct.p_a(idxA,:)) .* xStruct.c_H2_a(idxA,:);
    yStruct.n_dot_H2O_a_out = -2*p.K_a/p.delta_z*(u.p_a_out - xStruct.p_a(idxA,:)) .* xStruct.c_H2O_a(idxA,:);
    yStruct.a_H2O_a_out = (p.R * xStruct.T_a(idxA,:) .* xStruct.c_H2O_a(idxA,:)) ./ ...
    p.p_sat(xStruct.T_a(idxA,:));

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


% solid temperature/Humidity

yStruct.T_S = xStruct.T_s;
yStruct.Lambda = xStruct.lambda_m;
%average humidity
yStruct.a_H2O_avg = 0.5 * (yStruct.a_H2O_a_out + yStruct.a_H2O_c_out);


y_meas = [ ...
    yStruct.U_cell;
    yStruct.T_a_out;
    yStruct.T_c_out;
    yStruct.T_S;
    yStruct.p_a_out;
    yStruct.p_c_out;
    yStruct.a_H2O_a_out;
    yStruct.a_H2O_c_out ];

y_ctrl = [yStruct.T_S;
    yStruct.Lambda;
    yStruct.a_H2O_avg ];

end

