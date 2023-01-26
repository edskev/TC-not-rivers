function beta = plsq_beta_make(a,b,t)
%Makes the list of parameters beta
%   a       parameters for x equation
%   b       parameters for y eqution
%   t       t coordinates of the points
    
    a   = reshape(a,[],1);
    b   = reshape(b,[],1);
    t   = reshape(t,[],1);
    
    beta = [a;b;t];
end