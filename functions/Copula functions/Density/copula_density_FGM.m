function cuv=copula_density_FGM(u,v,a)

% parameter -1<=a<=1

cuv=1+a*(1-2*u).*(1-2*v);