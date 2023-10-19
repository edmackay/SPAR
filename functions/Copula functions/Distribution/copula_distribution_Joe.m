function C=copula_distribution_Joe(u,v,a)

% parameter a>=1
u1=(1-u).^a;
v1=(1-v).^a;
z = u1 + v1 - u1.*v1;
C = 1-z.^(1/a);