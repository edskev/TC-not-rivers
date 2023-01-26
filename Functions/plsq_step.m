function D_beta = plsq_step(beta,N,M,x,y,idx)
%Calculates the next least squares step
%   beta        Array of parameters [a;b;t]
%   N,M         Length of a and b
%   x,y         Data to fit to
%   idx         The indicies of beta that can change
%
%   D_beta      The recommended change in beta
    
    [rx,ry,F,G] = plsq_residual(beta,N,M,x,y);
    F  = F(:,idx);
    G  = G(:,idx);
    
    %Solve for D_beta
    D_beta  = zeros(size(beta));
    D_beta(idx)  = (F.'*F+G.'*G)\(F.'*rx+G.'*ry);
end