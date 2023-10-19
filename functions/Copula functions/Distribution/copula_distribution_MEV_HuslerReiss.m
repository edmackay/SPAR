function C=copula_distribution_MEV_HuslerReiss(u,v,a)

% parameter a>0
t=-log(u.*v);
w=-log(v)./t;

% z=log((1-w)./w); %% less numerically accurate
z=log(-log(u)) - log(-log(v));
x1=1/a + a/2.*z;
x2=1/a - a/2.*z;

term1=normcdf(x1);
term2=normcdf(x2);

A=(1-w).*term1 + w.*term2;

C=exp(-t.*A);

% x=-log(u);
% y=-log(v);
% t1=normcdf(1/a + a/2.*log(x./y));
% t2=normcdf(1/a + a/2.*log(y./x));
% 
% C=exp(-x.*t1 - y.*t2);