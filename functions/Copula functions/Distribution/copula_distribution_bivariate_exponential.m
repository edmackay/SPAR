function Cuv=copula_distribution_bivariate_exponential(u,v,a)

% parameter 0<=a<=1

Lu=log(u);
Lv=log(v);
Cuv = u.*v.* exp(-a *Lu.*Lv);