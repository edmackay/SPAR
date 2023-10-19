function C=copula_distribution_Nelsen16(u,v,a)

% parameter a>=0
z = u + v - 1 - a*(1./u + 1./v -1);
C = 0.5*(z+sqrt(z.^2 + 4*a));