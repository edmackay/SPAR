function C=copula_distribution_FGM(u,v,a)

% parameter a must be in [-1,1]

C=u.*v+a*u.*v.*(1-u).*(1-v);