function fx = margin_pdf(x,MARGIN)

if strcmpi(MARGIN.name,'logistic')
    fx = margin_pdf_logistic(x);
elseif strcmpi(MARGIN.name,'weibull')
    fx = margin_pdf_weibull(x,MARGIN.k);
elseif strcmpi(MARGIN.name,'laplace')
    fx = margin_pdf_laplace(x);
elseif strcmpi(MARGIN.name,'cauchy')
    fx = margin_pdf_cauchy(x);
elseif strcmpi(MARGIN.name,'Gaussian')    
    fx = normpdf(x);       
elseif strcmpi(MARGIN.name,'symGP')
    fx = margin_pdf_symGP(x,MARGIN.xi);
elseif strcmpi(MARGIN.name,'GP')
    fx = gppdf(x,MARGIN.xi,1,0);
elseif strcmpi(MARGIN.name(1:3),'exp')
    fx = exp(-x);
elseif strcmpi(MARGIN.name,'Pareto')
    fx = x.^(-2);
end
