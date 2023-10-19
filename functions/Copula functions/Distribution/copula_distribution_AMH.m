function C=copula_distribution_AMH(u,v,a)

% parameter -1<=a<=1

C = u.*v./(1-a*(1-u).*(1-v));