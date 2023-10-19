function x = margin_icdf(p,MARGIN)

if strcmpi(MARGIN.name,'logistic')
    x = margin_icdf_logistic(p);
elseif strcmpi(MARGIN.name,'laplace')
    x = margin_icdf_laplace(p);
elseif strcmpi(MARGIN.name,'cauchy')    
    x = margin_icdf_cauchy(p);
elseif strcmpi(MARGIN.name,'symGP')
    x = margin_icdf_symGP(p,MARGIN.xi);
end
