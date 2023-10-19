function C=copula_distribution_Clayton(u,v,a)

% parameter a>0

z=u.^(-a)+v.^(-a)-1;
C=z.^(-1/a);