function c=copula_density_Frank(u,v,a)

% parameter a in (-inf,inf), a~=0

g1=exp(-a)-1;
gu=exp(-a*u)-1;
gv=exp(-a*v)-1;
guv=exp(-a*(u+v))-1;
c=-a*g1*(guv+1)./(gu.*gv+g1).^2;