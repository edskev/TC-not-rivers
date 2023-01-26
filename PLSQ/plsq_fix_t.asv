function beta = plsq_fix_t(beta,N,M,x,y,idx)
%Finds a better value of t by analysing all values of t
%   beta        Array of parameters [a;b;t]
%   N,M         Length of a and b
%   x,y         Data to fit to
%   idx         Optional, the indicies in t that can change

    [a,b,t] = plsq_beta_split(beta,N,M);
    
    if nargin<6
        idx = 1:numel(t);
    else
        idx = reshape(idx,1,[]);
    end

    t_lin   = linspace(min(t),max(t),100);
    for i=idx
        rx0     = x(i) - plsq_poly(a,t(i));
        ry0     = y(i) - plsq_poly(b,t(i));
        
        rx      = x(i) - plsq_poly(a,t_lin);
        ry      = y(i) - plsq_poly(b,t_lin);
        
        r0      = rx0.^2 + ry0.^2;
        r       = rx.^2 + ry.^2;
        
        [~,rmin]= min(r);
        if r(rmin)<r0
            t(i)    = t_lin(rmin);
        end
    end
    
    beta = plsq_beta_make(a,b,t);
end