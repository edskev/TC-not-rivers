function [m,c]=plsq_tangent(beta,N,M,t0)
%Coefficients of the line y=mx+c passing through the curve at point t
%   beta        Array of parameters [a;b;t]
%   N,M         Length of a and b
%   t0          Point on curve to find tangent of
    
    [a,b] = plsq_beta_split(beta,N,M);
    
    x    = plsq_poly(a,t0);
    y    = plsq_poly(b,t0);
    dxdt = plsq_poly_diff_t(a,t0);
    dydt = plsq_poly_diff_t(b,t0);
    
    m = dydt./dxdt;
    c = y-m.*x;
end