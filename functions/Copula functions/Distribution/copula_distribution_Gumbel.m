function C=copula_distribution_Gumbel(u,v,a)

% parameter 0<a<=1

C = u.*v.*exp(-a*log(u).*log(v));