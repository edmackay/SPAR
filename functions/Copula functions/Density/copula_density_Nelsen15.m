function c=copula_density_Nelsen15(u,v,a)

% parameter a>=1
x=1-u.^(1/a);
y=1-v.^(1/a);
z = 1 - (x.^a + y.^a).^(1/a);

t1 = (u.*v).^(1/a - 1);
t2 = (x.*y).^(a - 1);
t3 = (x.^a + y.^a).^(1/a - 2);
c=(1-1/a) * t1 .* t2 .* t3 .* z.^(a - 2);
c(z<0)=0;