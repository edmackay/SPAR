function C=copula_distribution_Nelsen15(u,v,a)

% parameter a>=1
x=1-u.^(1/a);
y=1-v.^(1/a);
z = 1 - (x.^a + y.^a).^(1/a);
z(z<0)=0;
C = z.^a;