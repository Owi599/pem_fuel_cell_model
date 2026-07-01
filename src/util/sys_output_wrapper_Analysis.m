function y = sys_output_wrapper_Analysis(x,u,p)
    x = x(:);
    u = u(:);
    tmp.data = u.';
    [~,~,~,y] = sys_output_PEMFC(x,tmp,p);
    y = y(:);
end