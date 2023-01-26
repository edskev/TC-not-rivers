function f = plsq_poly(a,t)
%Computes f(a,t) or g(b,t) (using symmetry we just write code for f)
%   a           Array of ALL coefficients of the polynomial
%   t           Collumn vector of points
%   x           Collumn vector of polynmoial values
    
    t   = reshape(t,[],1);
    a   = reshape(a,1,[]);
    p   = 0:(numel(a)-1);
    
    f   = sum(  a.*t.^p  ,2);
end