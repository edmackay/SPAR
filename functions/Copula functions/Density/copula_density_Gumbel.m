function c=copula_density_Gumbel(u,v,a)

% parameter 0<a<=1

c = exp(-a*log(u).*log(v)).*(a^2 * log(u).*log(v) - a*(log(u)+log(v))-a+1);