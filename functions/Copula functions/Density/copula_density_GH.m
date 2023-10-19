function cuv=copula_density_GH(u,v,a)

% parameter 0<a<=1

cuv=(1 + a*(u + v + u.*v - 2) - a^2 *(u + v - u.*v - 1)) ./ (1 - a*(1-u).*(1-v)).^3;