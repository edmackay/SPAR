function C=copula_distribution_BB2(u,v,a,b)

% parameters a>0, b>0

C = (1 + (1/a) * log(exp(a*(u.^(-b)-1)) + exp(a*(v.^(-b)-1)) - 1)).^(-1/b);