function c=copula_density_BB9(u,v,params)

% parameters a>=1, b>=0
a=params(1);
b=params(2);

x=b-log(u);
y=b-log(v);
z=x.^a + y.^a - b^a;

C=exp(b-z.^(1/a));

term1=C./(u.*v);
term2=(x.*y).^(a-1);
term3=z.^(1/a-2);
term4=z.^(1/a) + a - 1;

c=term1.*term2.*term3.*term4;