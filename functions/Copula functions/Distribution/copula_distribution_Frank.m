function C=copula_distribution_Frank(u,v,a)

% parameter a in (-inf,inf), a~=0

g1=exp(-a)-1;
gu=exp(-a*u)-1;
gv=exp(-a*v)-1;
C=-(1/a)*log(1+gu.*gv/g1);