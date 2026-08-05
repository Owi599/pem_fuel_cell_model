function [mpcobj, mv_bounds_norm, range_u_norm] = design_mpc(sys, idxMV, idxMD, idxUD, mv_bounds, mu_u, sigma_u, Ts)
nu = numel(idxMV);
mv_bounds_norm = (mv_bounds - mu_u(idxMV)') ./ sigma_u(idxMV)';
range_u_norm = diff(mv_bounds_norm, 1, 2);

plant = setmpcsignals(sys, 'MV', idxMV, 'MD', idxMD, 'UD', idxUD);
plant.InputGroup.Unmeasured = idxUD;

mpcobj = mpc(plant, Ts, 200, 10);
mpcobj.Weights.OutputVariables = [1 1];
mpcobj.Weights.ManipulatedVariables = 0.01 * ones(1, nu);
mpcobj.Weights.ManipulatedVariablesRate = 1 ./ range_u_norm(:)';

for i = 1:nu
    mpcobj.MV(i).ScaleFactor = 0.05 * range_u_norm(i);
    mpcobj.MV(i).Min = mv_bounds_norm(i,1);
    mpcobj.MV(i).Max = mv_bounds_norm(i,2);
    mpcobj.MV(i).RateMin = -0.05 * range_u_norm(i);
    mpcobj.MV(i).RateMax =  0.05 * range_u_norm(i);
end
mpcobj.OV(1).ScaleFactor = 3;
mpcobj.OV(2).ScaleFactor = 10;
setEstimator(mpcobj, 'default');   % keep default here; switch to custom only for the nonlinear loop
review(mpcobj);
end