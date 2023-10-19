function cuv=copula_density_AMH(u,v,a)

% parameter -1<=a<=1

if a==1
    cuv=2*u.*v ./ (u+v-u.*v).^3;
else
    cuv=(1 + a*(u + v + u.*v - 2) - a^2 *(u + v - u.*v - 1)) ./ (1 - a*(1-u).*(1-v)).^3;
end