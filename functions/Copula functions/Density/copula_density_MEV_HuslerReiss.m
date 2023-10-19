function c=copula_density_MEV_HuslerReiss(u,v,a)

% parameter a>0

t=-log(u.*v);
w=-log(v)./t;

% z=log((1-w)./w); %% less numerically accurate
z=log(-log(u)) - log(-log(v));
x1=1/a + a/2.*z;
x2=1/a - a/2.*z;

term1=normcdf(x1);
term2=normcdf(x2);
term3=normpdf(x1);

A=(1-w).*term1 + w.*term2;

T1=exp(-t.*(A-1));
T2=term1.*term2;
T3=a*term3./(2*w.*t);

c=T1.*(T2 + T3);

c(isnan(c))=0;
c(isinf(c))=0;
