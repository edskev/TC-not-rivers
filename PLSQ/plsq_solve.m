function beta = plsq_solve(beta,N,M,x,y,step_prop,idx)
%Solves the least squares problem for beta
%   beta        Array of parameters [a;b;t]
%   N,M         Length of a and b
%   x,y         Data to fit to
%   step_prop   The scaling factor for the step size
%   idx         Optional, the indicies of beta that can change
    
    S = Inf;    %Step size
    
    if nargin<7
        idx = 1:numel(beta);
    end
    
    iter = 0;
    while S>1e-4 && iter<100
        D_beta  = plsq_step(beta,N,M,x,y,idx);
        S       = norm(D_beta)/numel(D_beta);
        beta    = beta + step_prop*D_beta;
        
        iter = iter+1;
    end
    
    plsq_plot(beta,N,M,x,y);
    drawnow;
end