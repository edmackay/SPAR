function c=copula_density_Nelsen12(u,v,a)

% parameter a>=1

uinv=1./u;
vinv=1./v;
t1=(uinv - 1).^a;
t2=(vinv - 1).^a;
t3=(uinv - 1).^(a-1);
t4=(vinv - 1).^(a-1);
z=t1+t2;

term1=(uinv.*vinv).^2;
term2=(z.^(1/a)+ 1).^(-3);
term3=z.^(1/a - 2);
term4=(a+1)* z.^(1/a)+ a - 1;
c=term1.*t3.*t4.*term2.*term3.*term4;
