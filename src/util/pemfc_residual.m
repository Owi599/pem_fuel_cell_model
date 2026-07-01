function res = pemfc_residual(tt, x, xdot, u_fun, p, M)
    u = u_fun(tt);
    u = u(:);                % force column vector
    rhs = ode_PEMFC(tt, x, u);
    res = M*xdot - rhs(:);   
end