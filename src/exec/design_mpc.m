function [mpcobj, mv_bounds_norm, range_u_norm] = ...
    design_mpc(sys, idxMV, idxMD, idxUD, mv_bounds, mu_u, sigma_u, Ts)

nu = numel(idxMV);

mv_bounds_norm = (mv_bounds - mu_u(idxMV).') ./ sigma_u(idxMV).';
range_u_norm = mv_bounds_norm(:,2) - mv_bounds_norm(:,1);

plant = setmpcsignals(sys, ...
    'MV', idxMV, ...
    'MD', idxMD, ...
    'UD', idxUD);

mpcobj = mpc(plant, Ts,30 , [2 3 5 10 30 50]);

% Already normalized model: preserve unit scale factors.
mpcobj.OV(1).ScaleFactor = 1;
mpcobj.OV(2).ScaleFactor = 1;

for i = 1:nu
    mpcobj.MV(i).ScaleFactor = range_u_norm(i);

    % Hard physical actuator limits.
    mpcobj.MV(i).Min = mv_bounds_norm(i,1);
    mpcobj.MV(i).Max = mv_bounds_norm(i,2);

    % Initial actuator slew limits: ±1% of permitted range per sample.
    mpcobj.MV(i).RateMin = -0.005 * range_u_norm(i);
    mpcobj.MV(i).RateMax =  0.005 * range_u_norm(i);

    % Preserve physical magnitude limits as hard constraints.
    mpcobj.MV(i).MinECR = 0;
    mpcobj.MV(i).MaxECR = 0;

    % Let rate constraints relax if necessary to return within hard MV bounds.
    mpcobj.MV(i).RateMinECR = 1;
    mpcobj.MV(i).RateMaxECR = 1;
end

% Equal initial importance for temperature and humidity tracking.
mpcobj.Weights.OutputVariables = [1,1];

% Set to zero initially to avoid deliberate output offset.
mpcobj.Weights.ManipulatedVariables = zeros(1, nu);

% Main damping parameter.
mpcobj.Weights.ManipulatedVariablesRate = [2.5000 10 2.5000 5 0.2500];  % heavier damping on MV2, MV4

setEstimator(mpcobj, 'default');

review(mpcobj);

end