function c=copula_density_Joe(u,v,a)

% parameter a>=1
u1=(1-u).^a;
v1=(1-v).^a;
z = u1 + v1 - u1.*v1;
c = ((1-u).*(1-v)).^(a-1) .* z.^(1/a - 2) .* (z+a-1);
c(isnan(c))=0;

