function [rx,ry,F,G] = plsq_residual(beta,N,M,x,y)
%Calculates the residuals and derivatives for the least squares method
%   beta        Array of parameters [a;b;t]
%   N,M         Length of a and b
%   x,y         Data to fit to
    
    [a,b,t] = plsq_beta_split(beta,N,M);
    I       = numel(t);
    
    %Residuals
    rx      = x - plsq_poly(a,t);
    ry      = y - plsq_poly(b,t);
    
    %Derivatives of f and g wrt beta
    F       = [ plsq_poly_diff_a(a,t)   , zeros(I,M+1)          , plsq_poly_diff_t(a,t) ];
    G       = [ zeros(I,N+1)            , plsq_poly_diff_a(b,t) , plsq_poly_diff_t(b,t) ];
end