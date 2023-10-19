function C=copula_distribution_MEV_symlog(u,v,a)

% parameter a>=1

Lu=-log(u);
Lv=-log(v);
z=Lu.^a + Lv.^a;

C=exp(-z.^(1/a));