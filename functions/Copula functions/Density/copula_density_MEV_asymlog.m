function c=copula_density_MEV_asymlog(u,v,params)

% parameters
% a>=1
% 0<=psi1, psi2<=1

a=params(1);
psi1=params(2);
psi2=params(3);

t=-log(u)-log(v);
w=-log(v)./t;

z = (psi1*(1-w)).^a + (psi2*w).^a;
A = (1-psi1)*(1-w) + (1-psi2)*w + z.^(1/a);
A1 = psi1 - psi2 + (psi2^a * w.^(a-1) - psi1^a *(1-w).^(a-1)) .* z.^(1/a-1);
A2 = (a-1) * (psi1*psi2)^a * (w.*(1-w)).^(a-2) .* z.^(1/a-2);

T1=exp(-t.*(A-1));
T2=(A + (1-w).*A1).*(A - w.*A1);
T3=w.*(1-w)./t;

c=T1.*(T2 + T3.*A2);

c(isnan(c))=0;
c(isinf(c))=0;

