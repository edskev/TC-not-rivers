function [a,b,t] = plsq_beta_split(beta,N,M)
%Splits up the list of parameters beta
%   beta        Array of parameters [a;b;t]
%   N,M         Length of a and b
    
    beta = reshape(beta,[],1);
    
    a = beta(1:N+1).';
    b = beta(N+1+(1:M+1)).';
    t = beta(N+M+3:end);
end