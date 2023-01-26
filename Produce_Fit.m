%This script fits two functions, x(t) and y(t), to a data set of pairs (x,y).
%The functions are of the following forms:
%       x = f(a,t) = a_0 + a_1 t^1 + ... + a_(N-1) t^(N-1) + a_N t^N
%       y = g(a,t) = b_0 + b_1 t^1 + ... + b_(M-1) t^(M-1) + b_M t^M
%Note:  the restriction on coefficients is to remove gauge symmetry.
%Note:  the form of the expressions requires N>=1 and M>=1.

N=3;
M=2;

%----Load in the data for fluvial-----
data = csvread('OC_withoutTsunami.csv',1,1);
x = log10(data(:,2));
y = log10(data(:,1));

t               = ones(size(x))*-1;
t(x>-2 & y<0)   = 1;
t(y>0)          = 0;
t               = t + (rand(size(t))-0.5)*0.1;

% %----Load in the data for GC, log law production----
% data = csvread('GC_Pll.csv',1,1);
% x = log10(data(:,2));
% y = log10(data(:,1));
% 
% t               = ones(size(x))*-1;
% t(x>-2 & y<1)   = 1;
% t(y>1)          = 0;
% t               = t + (rand(size(t))-0.5)*0.1;

% %----Load in the data for GC, top hat production----
% data = csvread('GC_PTBT.csv',1,1);
% x = log10(data(:,2));
% y = log10(data(:,1));
% 
% t               = zeros(size(x));
% t(y<-2)         = -1;
% t(y<-3)         = -2;
% t(y>-1)         = 1;
% t               = t + (rand(size(t))-0.5)*0.1;


%Initialise a and b
a = zeros(1,N+1); a(2) = 1;
b = rand(1,M+1);  b(2) = 0;
beta = plsq_beta_make(a,b,t);

plsq_plot(beta,N,M,x,y);
drawnow;

%Solve for b
idx  = [1,N+2,N+4:N+M+2];
beta = plsq_solve(beta,N,M,x,y,1,idx);

%Solve for t
idx =  [N+M+3:numel(beta)];
beta = plsq_solve(beta,N,M,x,y,0.1,idx);

%Solve for a,b then fix t
idx =  [1,3:N+2,N+4:numel(beta)];
beta = plsq_solve(beta,N,M,x,y,0.1,idx);
for i=1:2
    beta = plsq_fix_t(beta,N,M,x,y);
    beta = plsq_solve(beta,N,M,x,y,0.1,idx);
end

%If dx/dt<0 anywhere then set min(dx/dt)=0 and solve
[a,b,t] = plsq_beta_split(beta,N,M);
if a(2)-a(3)^2/(3*a(4))<0
    idx =  [1,3,N+2,N+4:numel(beta)];
    
    for i=1:2
        beta = plsq_fix_t(beta,N,M,x,y,find(y>0));
        S = Inf;
        
        while S>1e-2 || (S>1e-4 && p==0)
            beta_old = beta;
            
            [rx,ry,F,G] = plsq_residual(beta,N,M,x,y);

            da4da       = [0 , 0 , 2*a(3)/3 , NaN];
            F(:,1:N+1)  = F(:,1:N+1) + F(:,4).*da4da;

            F  = F(:,idx);
            G  = G(:,idx);

            %Solve for D_beta
            D_beta          = zeros(size(beta));
            D_beta(idx)     = (F.'*F+G.'*G)\(F.'*rx+G.'*ry);
            beta            = beta + 0.1*D_beta;
            
            [a,b,t] = plsq_beta_split(beta,N,M);
            a(4) = (a(3)^2)/3;
            beta = plsq_beta_make(a,b,t);
            
            S = 10*norm(beta-beta_old)/numel(beta);

        end
        plsq_plot(beta,N,M,x,y);
        drawnow;
    end
end

a
b
T = [min(t),max(t)]