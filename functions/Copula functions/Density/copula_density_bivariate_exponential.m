function cuv=copula_density_bivariate_exponential(u,v,a)

% parameter 0<=a<=1

Lu=log(u);
Lv=log(v);
T1=exp(-a *Lu.*Lv);
T2=(a^2 *Lu.*Lv - a*(Lu + Lv) - a + 1);
cuv = T1 .*T2;