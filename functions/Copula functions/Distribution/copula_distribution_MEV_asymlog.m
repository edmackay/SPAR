function C=copula_distribution_MEV_asymlog(u,v,params)

% parameter a>=1

a=params(1);
psi1=params(2);
psi2=params(3);

Lu=-log(u);
Lv=-log(v);
z=(psi1*Lu).^a + (psi2*Lv).^a;

C=exp(-(1-psi1)*Lu - (1-psi2)*Lv - z.^(1/a));