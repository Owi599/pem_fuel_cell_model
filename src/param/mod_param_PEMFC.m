function out = mod_param_PEMFC()
%mod_param_PEMFC generates default model parameter structure
%               -no fields added to p after p is read for the first time
%               -do not re-create p if p already exists
%               -adaptable parameters need to be complex for complex step
%                method in generated code

%% declare persistent variables
% persistent p

%% get other parameter structs


%% initialize parameter values

% if isempty(p)

% --------------------------------------------------------------------------------------  
% stack geometry
p.N_cells = 20; % number of cells in the stack
p.A_cell = 284e-4; % active cell area in m^2

% --------------------------------------------------------------------------------------  
% geometry data of the individual layers 
L_z = 0.2462;              % channel length [m]
p.L_z=L_z;
p.L_x = 1.1e-3;            % channel height [m] (A_cell/(L_z*N_channels))
p.delta_a = 0.12e-3;       % anode channel width [m] (approx.)
p.delta_a_g = 0.15e-3;     % anode GDL width [m] 0.12-0.215 mm
p.delta_a_c = 0.01e-3;     % anode CL width [m]
p.delta_m = 0.012134e-3;   % membrane width [m] 0.012-0.018 mm
p.delta_s = 0.4e-3;        % solid width [m] approx. (sum CL + GDL + M)
p.delta_c_c = 0.01e-3;     % cathode CL width [m]
p.delta_c_g = 0.151e-3;    % cathode GDL width [m] 0.12-0.215 mm
p.delta_c = 0.172e-3;      % cathode channel width [m]  (approx.)

% --------------------------------------------------------------------------------------  
% material data of the membrane and the solid [Neubr99]
x_m_dry = 0.909;            % ion exchange capacity [mol kg^-1]
p.x_m_dry =x_m_dry;
rho_m_dry = 2.05e3;         % membrane density [kg m^-3]
p.rho_m_dry = rho_m_dry ;
rho_H2O_l = 1e3;            % water density [kg m^-3]
p.rho_H2O_l = rho_H2O_l;
p.lambda_s = 0.43;          % thermal conductivity [W m^-1 K^-1]
p.rho_c_p_m = 2.31e7;        % [J m^-3 K^-1]  
p.rho_c_p_s = 4.15e6;        % [J m^-3 K^-1]

% --------------------------------------------------------------------------------------  
% material data of the reactants and products [Perry97, VDI06]
p.D_eff_H2_H2O_g = 1e-6;    % diffusion coefficient of the pair H2 and H2O (g) [m^2 s^-1]
p.D_eff_O2_H2O_g = 2.55e-6;    % diffusion coefficient of the pair O2 and H2O (g) [m^2 s^-1]
p.D_eff_O2_N2 = 1.13e-6;       % diffusion coefficient of the pair O2 and N2 [m^2 s^-1]
p.D_eff_N2_H2O_g = 2.64e-6;  % diffusion coefficient of the pair N2 and H2O (g) [m^2 s^-1]
p.mu_0_H2O_g = -228.597;    % standard electrochemical potential of H2O (g) [J mol^-1]
p.c_p_H2 = 28.824;          % specific heat capacity of H2 [J (mol*K)^-1]
p.c_p_O2 = 29.335;          % specific heat capacity of O2 [J (mol*K)^-1]
p.c_p_N2 = 29.125;          % specific heat capacity of N2 [J (mol*K)^-1]
p.c_p_H2O_g = 33.58;        % specific heat capacity of H2O (g) [J (mol*K)^-1]
p.c_p_H2O_l = 75.312;       % specific heat capacity of H2O (l) [J (mol*K)^-1]
p.h_H2_ref = 0;             % molar enthalpy [J mol^-1] for reference temperature 298.15 [K]
p.h_H2O_g_ref = -241830;    % molar enthalpy [J mol^-1] for reference temperature 298.15 [K]
p.h_O2_ref = 0;             % molar enthalpy [J mol^-1] for reference temperature 298.15 [K]
p.h_N2_ref = 0;             % molar enthalpy [J mol^-1] for reference temperature 298.15 [K]
p.M_O2 = 32e-3;             % oxygen molar mass [kg/mol]
p.M_N2 = 28e-3;             % nitrogen molar mass [kg/mol]
M_H2O = 18e-3;              % water molar mass [kg/mol]
p.M_H2O=M_H2O;
p.M_H2 = 2.016e-3;          % hydrogen molar mass [kg/mol]

% --------------------------------------------------------------------------------------  
% further system parameters [Baern92, Peng08, VDI06, Woehr99]
p.R = 8.31446261815324;     % universal gas const. [J (mol*K)^-1]
F = 9.648533212331e4;       % Faraday constant [C mol^-1]
p.F = F;
T_ref = 298.15;             % reference temperature for gases [K]
p.T_ref =T_ref;
T_ref_s = 353.15;           % reference temperature for solid [K]
p.T_ref_s =T_ref_s;
p.p_u = 101325;             % atm. pressure [Pa]
p.alpha_1 = 96.5;           % heat transfer coefficient between solid and gas channels [W m^-2 K^-1]
p.alpha_2 = 1524;           % heat transfer coefficient between coolant and solid [W m^-2 K^-1] % !!! Mangold2010 100
p.lambda_a = 0.189;         % heat conduction coefficient of anode [W m^-1 K^-1]
p.lambda_c = 0.288;         % heat conduction coefficient of cathode [W m^-1 K^-1]

% --------------------------------------------------------------------------------------  
% enthalpy calculation: f(enthalpy for reference temperature, specific heat capacity, temperature) 
p.h_T =  @(h_x_ref, c_p_x, T) h_x_ref + c_p_x .* (T - T_ref);

% --------------------------------------------------------------------------------------  
p.f_v = 60;                 % surface correction factor
p.i_0_ref_a = 104;          % anode reference current density [A m^-2]
p.i_0_ref_c = 1e-5;         % cathode reference current density [A m^-2]
p.C_a = 8.25e6;             % anode double layer capacity [F m^-3]
p.C_c = 8.25e6;             % cathode double layer capacity [F m^-3]
p.Delta_g_0 = 73000;        % activation energy of cathode reaction [J mol^-1]
p.K_a = 7.69e-5;            % anode Ergun coefficient [m^2 s^-1 Pa^-1]
p.K_c = 9.33e-5;            % cathode Ergun coefficient [m^2 s^-1 Pa^-1]
p.p_H2_ref = 101325;        % reference pressure (anode) [Pa]
p.p_O2_ref = 101325;        % reference pressure (cathode) [Pa]
p.alpha_a = 0.95;           % transfer coefficient (anode)
p.alpha_c = 0.5;            % transfer coefficient (cathode)
delta_G = -228.6e3;         % Gibbs free energy (H2O(gas)) at 25 C, 1 atm
p.delta_G = delta_G;
p.E_cr = 66e3;              % molar activation energy of the cathode reaction [J/mol]

% --------------------------------------------------------------------------------------  
% anode open circuit potential (assumed to be zero)
p.delta_Phi_a_ref = 0;

% cathode open circuit potential (Nernst equation)
p.delta_Phi_c_ref = -delta_G/(2*F);

% additional gains
p.k_lambda = 1;
p.k_D_w = 1;
p.k_kappa = 1;
p.k_t = 1;

% --------------------------------------------------------------------------------------  
% membrane model parameters [Neubr99] 
% transport number of water: t_w = f(water content)
k_t_w = 0.9222 * 2.5/22;
p.k_t_w = k_t_w;
p.t_w = @(lambda_m) k_t_w * lambda_m;

% diffusion coefficient of water: D_w = f(water content, temperature) 
T_ref_D_w = 353.15;
if 0
    % initial model
    p.D_w = @(lambda_m, T_s) 10.^(-10.775 + 0.3436 * lambda_m - 0.0189 * lambda_m.^2 + 0.0004 * lambda_m.^3)...
            .* exp(-(2640 * exp(-0.6 * lambda_m) + 1517) .* (1 ./ T_s - 1 / T_ref_D_w));
else
    % simplified model - linear appox.
    k_D1 = 2.24574e-10;
    k_D2 = 1.796592e-10; % if k_D1=k_D2 -> no temperature dependency
    p.D_w=@(lambda_m,T_s)  k_D1*lambda_m.*(k_D2/k_D1 + (T_s-273.15-50)/(80-50) * (1-k_D2/k_D1));
end



% conductivity of the membrane: kappa = f(water content, temperature)
T_ref_kappa = 353.15;
if 0
    % initial model
    p.kappa = @(lambda_m, T_s) (0.0013 * lambda_m.^3 + 0.0298 * lambda_m.^2 + 0.2658 * lambda_m) .* ...
            exp(-(2640 * exp(-0.6 * lambda_m) + 1183) .* (1./T_s - 1 / T_ref_kappa));
else
    % simplified model - linear approx.
    k_k1 = 0.45; % 80 °C
    p.kappa=@(lambda_m, T_s) k_k1 * lambda_m;
end
    
% membrane density [kg m^-3]: rho_m = f(water content)
p.rho_m = @(lambda_m) rho_H2O_l * rho_m_dry * ...
    (1 + lambda_m * M_H2O * x_m_dry)./(rho_H2O_l + lambda_m * M_H2O * x_m_dry * rho_m_dry);

% ion exchange capacity [mol/kg]: x_m = f(water content)
p.x_m = @(lambda_m) 1./(1 + lambda_m * M_H2O * x_m_dry) * x_m_dry;

% mole fraction of water in membrane: xi_H2O = f(water content)
p.xi_H2O = @(lambda_m) lambda_m ./ (lambda_m + 2);

% mole fraction of protones in membrane: xi_H = f(water content)
p.xi_H = @(lambda_m) 1 ./ (lambda_m + 2);

% saturation pressure of H20 [Pa]: p_sat = f(vapor_temperature) 
p.p_sat = @(T) (100 * exp(19.016 - 4064.95 ./ (T - 273.15 + 236.25))); 

% water content on the boundaries of membrane: lambda_m_x = f(relative humidity, temperature)
% linear interpolated for temperatures t = 30 ... 80 C
k_lambda1 = 200; % correction factor (1-exp(-k_lambda*a_H2O)) for lambda(a=0,T)=0
p.lambda_m_x = @(a_H2O,T_s) ( (36 * a_H2O.^3 - 39.85 * a_H2O.^2 + 17.81 * a_H2O + 0.043) ...
                             +((T_s - 273.15) - 30).*( (14.1 * a_H2O.^3 - 16 * a_H2O.^2 + 10.8 * a_H2O + 0.3)...
                                                      -(36 * a_H2O.^3 - 39.85 * a_H2O.^2 + 17.81 * a_H2O + 0.043) ...
                                                     ) ./ (80 - 30)  ...
                            ) .* (1-exp(-k_lambda1*a_H2O));

% --------------------------------------------------------------------------------------  
% discretization along the z-coordinate (N = 3 ... 50)
N = 20;               % number of discretization points
p.N=N;
delta_z = L_z / N;    % discretization step [m]
p.delta_z = delta_z;
% vector of spatial coordinate
zv=((1:N)-0.5)*delta_z;
p.zv=zv;

% --------------------------------------------------------------------------------------
% counter- or co-flow feeding in the gas channels (conter-flows used in experiments)
% direction is changed only in the anode channel 
p.counter_flow = 1;     % 1 - counter-flow feeding       
                        % 0 - co-flow feeding

% --------------------------------------------------------------------------------------  
% define states and initialize anonymous converter functions    
% state description: name, dimensions, and dynamic/algebraic
states.variableNames =  {     'c_H2_a';
                              'c_H2O_a';
                              'c_O2_c';
                              'c_H2O_c';...
                              'c_N2_c';...
                              'lambda_m';...
                              'delta_Phi_a';...
                              'delta_Phi_c';...
                              'rho_u_a';...
                              'rho_u_c';...
                              'rho_e_s';...
                              'xi_H2_ca';...
                              'xi_H2O_ca';...
                              'xi_O2_cc';...
                              'xi_N2_cc';...
                              'xi_H2O_cc';...
                              'delta_Phi_m';...
                              'i_m';...
                              'T_a';...
                              'T_c';...
                              'T_s';...
                              'p_a';...
                              'p_c';...
                              'U_cell'};
states.variableDimensions = [N;N;N;N;N;N;N;N;N;N;N;N;N;N;N;N;N;N;N;N;N;N;N;1];
states.dynamicState =       [1;1;1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0];

p.states=states;

p.state2struct = @(xvec)    vec2struct(xvec   ,states.variableNames,states.variableDimensions);
p.struct2state = @(xstruct) struct2vec(xstruct,states.variableNames,states.variableDimensions);

ndynamic = sum(states.variableDimensions(states.dynamicState==1));
nstates = sum(states.variableDimensions);
p.nstates = nstates;
p.ndynamic = ndynamic;

if ~coder.target('MATLAB')
    rhs_init = vec2struct((complex(zeros([nstates,1]))),states.variableNames,states.variableDimensions);% CK: complex for complex step method in generated code
    p.rhs_init = rhs_init;
end

% --------------------------------------------------------------------------------------
% define test bench input variables
% state description: name, dimensions
testbench.variableNames = {'p_So_C';...
                           'FN_Si_Air_C';...
                           'DPT_Si_C';...
                           'T_Si_C';...
                           'p_So_A';...
                           'FN_Si_H2_A';...
                           'DPT_Si_A';...
                           'T_Si_A';...
                           'T_Si_CL';...
                           'FN_Si_CL';...
                           'I_S';...
                           'p_Si_A';...
                           'p_Si_C'};
testbench.variableDimensions = [1;1;1;1;1;1;1;1;1;1;1;1;1];

p.testbench = testbench;

p.testbench2struct = @(xvec) vec2struct(xvec, testbench.variableNames, testbench.variableDimensions);

% --------------------------------------------------------------------------------------
% define model input variables
% state description: name, dimensions

inputs.variableNames = {'n_dot_H2_a_in';...
                        'n_dot_H2O_a_in';...
                        'n_dot_O2_c_in';...
                        'n_dot_H2O_c_in';...
                        'n_dot_N2_c_in';...
                        'T_a_in';...
                        'T_c_in';...
                        'T_cool';...
                        'p_a_out';...
                        'p_c_out';...
                        'I_cell'};
inputs.variableDimensions = [1;1;1;1;1;1;1;1;1;1;1];

p.inputs = inputs;

p.inputs2struct = @(xvec)    vec2struct(xvec,    inputs.variableNames,inputs.variableDimensions);
p.inputs2vec    = @(xstruct) struct2vec(xstruct, inputs.variableNames,inputs.variableDimensions);

% --------------------------------------------------------------------------------------  
% define fluxes between compartments and intialization struct
fluxes.variableNames = {'n_dot_H2O_am';...
                        'n_dot_H2O_cm';...
                        'n_dot_H2_a';...
                        'n_dot_H2O_a';...
                        'n_dot_O2_c';...
                        'n_dot_H2O_c';...
                        'n_dot_N2_c'};
fluxes.variableDimensions = [N;N;N;N;N;N;N];
p.fluxes = fluxes;

nfluxes = sum(fluxes.variableDimensions);

if ~coder.target('MATLAB')
    fluxes_init = vec2struct((complex(zeros([nfluxes,1]))),fluxes.variableNames,fluxes.variableDimensions);% CK: complex for complex step method in generated code
    p.fluxes_init = fluxes_init;
end

% --------------------------------------------------------------------------------------  
% mass matrix - constructed from information in states, diagonal entries
% only, states have to be sorted by dynamic first and algebraic last
% for DAE system xt = f(t,x), 0 = g(t,x) <=> M*xt=rhs(t,x) with M = [I 0; 0 0]

p.M = @() sparse(1:ndynamic,1:ndynamic,ones([ndynamic,1]),nstates,nstates);

% --------------------------------------------------------------------------------------
% define model outputs (not test bench outputs!)
% state description: name, dimensions
outputs.variableNames = {'U_cell';...
                        'T_a_out';...
                        'T_c_out';...
                        'n_dot_H2_a_out';...
                        'n_dot_H2O_a_out';...
                        'n_dot_O2_c_out';...
                        'n_dot_H2O_c_out';...
                        'p_a_out';...
                        'p_c_out';...
                        'a_H2O_a_out';...
                        'a_H2O_c_out'};
outputs.variableDimensions = [1;1;1;1;1;1;1;1;1;1;1];

p.outputs = outputs;

p.struct2output = @(xstruct) struct2vec(xstruct,outputs.variableNames,outputs.variableDimensions);
p.output2struct = @(xvec) vec2struct(xvec,outputs.variableNames,outputs.variableDimensions);

% output scaling - min-max-scaling
output2scaled = @(x, x_min, x_max) rescaling(x, x_min, x_max, 0, 1);
p.output2scaled = output2scaled;

% -------------------------------------------------------------------------------------- 
% selected outputs for the optimization of the hybrid model
p.selectedOutputs = [1 0 0 0 1 0 1 0 0 0 0]; % 1 - used | 0 - unused output

% --------------------------------------------------------------------------------------
% test bench outputs for identification
% state description: name, dimensions
tb_outputs.variableNames = {'U_S';...
                        'T_So_A';...
                        'T_So_C';...
                        'p_Si_A';...
                        'p_Si_C'};
tb_outputs.variableDimensions = [1;1;1;1;1];

p.tb_outputs = tb_outputs;

p.tb_struct2output = @(xstruct) struct2vec(xstruct,tb_outputs.variableNames,tb_outputs.variableDimensions);
p.tb_output2struct = @(xvec) vec2struct(xvec,tb_outputs.variableNames,tb_outputs.variableDimensions);

% output scaling - min-max-scaling
tb_output2scaled = @(x, x_min, x_max) rescaling(x, x_min, x_max, 0, 1);
p.tb_output2scaled = tb_output2scaled;

% -------------------------------------------------------------------------------------- 
% selected outputs for the optimization of the hybrid model
p.tb_selectedOutputs = [1 0 0 1 1]; % 1 - used | 0 - unused output


% end


out = p;

