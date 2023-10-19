function c=copula_density_MEV_symlog(u,v,a)

% parameter a>=1

tau=-log(u.*v);
psi=-log(v)./tau;

z = psi.^a + (1-psi).^a;
A = z.^(1/a);
T1 = exp(-tau.*(A-1));
T2 = z.^(2/a-2) .* (psi.*(1-psi)).^(a-1);
T3 = (a-1)*z.^(-1/a).*T2./tau;
c=T1.*(T2 + T3);

c(isnan(c))=0;
c(isinf(c))=0;

