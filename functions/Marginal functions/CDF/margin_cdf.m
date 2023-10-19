function Fx = margin_cdf(x,MARGIN)

if strcmpi(MARGIN.name,'logistic')
    Fx = margin_cdf_logistic(x);
elseif strcmpi(MARGIN.name,'laplace')
    Fx = margin_cdf_laplace(x);
elseif strcmpi(MARGIN.name,'cauchy')    
    Fx = margin_cdf_cauchy(x);
elseif strcmpi(MARGIN.name,'Gaussian')    
    Fx = normcdf(x);    
elseif strcmpi(MARGIN.name,'GP')
    Fx = gpcdf(x,MARGIN.xi,1,0);
elseif strcmpi(MARGIN.name,'symGP')
    Fx = margin_cdf_symGP(x,MARGIN.xi);
elseif strcmpi(MARGIN.name,'weibull')
    Fx = margin_cdf_weibull(x,MARGIN.k);    
elseif strcmpi(MARGIN.name(1:3),'exp')
    Fx = 1-exp(-x);    
elseif strcmpi(MARGIN.name,'Pareto')
    Fx = 1-x.^(-1);
end
