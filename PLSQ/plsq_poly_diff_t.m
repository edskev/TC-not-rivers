function dfdt = plsq_poly_diff_t(a,t)
%Computes df/da(a,t) or dg/db(b,t) (using symmetry we just write code for f)
%   a           Array of ALL coefficients of the polynomial
%   t           Collumn vector of points
%   x           Collumn vector of polynmoial values
    
    t   = reshape(t,[],1);
    a   = reshape(a,1,[]);
    p   = 0:(numel(a)-1);
    
    dfdt= diag(sum(  a(2:end).*p(2:end).*t.^p(1:end-1)  ,2));
end