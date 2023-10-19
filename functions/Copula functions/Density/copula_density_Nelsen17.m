function c=copula_density_Nelsen17(u,v,a)

% parameter -inf<a<inf, a~=0

t1=(1+u).^(-a) - 1;
t2=(1+v).^(-a) - 1;
S=t1.*t2;
z=1+S/(2^(-a)-1);

term1=((u+1).*(v+1)).^(-a-1);
term2=z.^(-1/a-1);
term3=(1/a + 1)*(S./z) + 1 -2^(-a);

c=(a/(2^(-a)-1)^2) .*term1.*term2.*term3;
