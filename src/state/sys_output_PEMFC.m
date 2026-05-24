function y = sys_output_PEMFC(x,u_traj_pp,p)
% --------------------------------------------------------------------------------------   
% model output function y = g(t,x,u)
% for assumed output vector (not test bench outputs!)
% y = [U_cell, T_a_out, T_c_out; n_dot_H2_a_out, n_dot_H2O_a_out,
%      n_dot_O2_c_out, n_dot_H2O_c_out, p_a_out, p_c_out];
% fnc inputs
    % x         - state vector
    % u_traj_pp - model input trajectories
    % p         - parameter structure
% fnc outputs
    % y         - output vector


% -------------------------------------------------------------------------
% output calculation from state vector and inputs
xStruct = p.state2struct(x);
u = p.inputs2struct(u_traj_pp.data.');

% voltage
yStruct.U_cell = xStruct.U_cell;
% output flows, temperature and pressure in the anode channel
if p.counter_flow
    % counter-flow case
    yStruct.n_dot_H2_a_out = -2*p.K_a/p.delta_z*(u.p_a_out - xStruct.p_a(1,:)) .* xStruct.c_H2_a(1,:);
    yStruct.n_dot_H2O_a_out = -2*p.K_a/p.delta_z*(u.p_a_out - xStruct.p_a(1,:)) .* xStruct.c_H2O_a(1,:);
    yStruct.T_a_out = xStruct.T_a(1,:);
    yStruct.p_a_out = xStruct.p_a(p.N,:);
    yStruct.a_H2O_a_out = (p.R * xStruct.T_a(1,:) .* xStruct.c_H2O_a(1,:)) ./ ...
            p.p_sat(xStruct.T_a(1,:));
else
    % co-flow case
    yStruct.n_dot_H2_a_out = -2*p.K_a/p.delta_z*(u.p_a_out - xStruct.p_a(p.N,:)) .* xStruct.c_H2_a(p.N,:);
    yStruct.n_dot_H2O_a_out = -2*p.K_a/p.delta_z*(u.p_a_out - xStruct.p_a(p.N,:)) .* xStruct.c_H2O_a(p.N,:);
    yStruct.T_a_out = xStruct.T_a(p.N,:);
    yStruct.p_a_out = xStruct.p_a(1,:);
    yStruct.a_H2O_a_out = (p.R * xStruct.T_a(p.N,:) .* xStruct.c_H2O_a(p.N,:)) ./ ...
        p.p_sat(xStruct.T_a(p.N,:));
end
% output flows in the cathode channel
yStruct.n_dot_O2_c_out = -2*p.K_c/p.delta_z*(u.p_c_out - xStruct.p_c(p.N,:)) .* xStruct.c_O2_c(p.N,:);
yStruct.n_dot_H2O_c_out = -2*p.K_c/p.delta_z*(u.p_c_out - xStruct.p_c(p.N,:)) .* xStruct.c_H2O_c(p.N,:);
yStruct.n_dot_N2_c_out = -2*p.K_c/p.delta_z*(u.p_c_out - xStruct.p_c(p.N,:)) .* xStruct.c_N2_c(p.N,:); % not used in the output vector
% temperature and pressure in the cathode channel
yStruct.T_c_out = xStruct.T_c(p.N,:);
yStruct.p_c_out = xStruct.p_c(1,:);
yStruct.a_H2O_c_out = (p.R * xStruct.T_c(p.N,:) .* xStruct.c_H2O_c(p.N,:)) ./ ...
    p.p_sat(xStruct.T_c(p.N,:));


% model output
y = p.struct2output(yStruct).';

end

