function C=copula_distribution_Nelsen17(u,v,a)

% parameter -inf<a<inf, a~=0

t1=(1+u).^(-a) - 1;
t2=(1+v).^(-a) - 1;
bot=2^(-a) - 1;
C = (1+t1.*t2/bot).^(-1/a) - 1;