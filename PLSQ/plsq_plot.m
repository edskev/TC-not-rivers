function plsq_plot(beta,N,M,x,y)
%Plots found curve from beta
%   beta        Array of parameters [a;b;t]
%   N,M         Length of a and b
%   x,y         Data to fit to
    
    [a,b,t] = plsq_beta_split(beta,N,M);
    
    curve_t = linspace(min(t),max(t),100);
    curve_x = plsq_poly(a,curve_t);
    curve_y = plsq_poly(b,curve_t);
    
    clf;
    hold on;
    plot(x,y,'xb','LineWidth',2)
    
    if nargin>3
        for i=1:numel(x)
            plot([x(i),plsq_poly(a,t(i))],[y(i),plsq_poly(b,t(i))],'-','Color',[0.7,0,0.7])
        end
    end

    plot(curve_x,curve_y,'-','Color',[0,0.7,0.7],'LineWidth',2);
    
    daspect([1 1 1])
end