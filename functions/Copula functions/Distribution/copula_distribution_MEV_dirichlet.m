function C=copula_distribution_MEV_dirichlet(u,v,params)

% parameter 0<a<=1

a=params(1);
b=params(2);

Lu=-log(u);
Lv=-log(v);
r=Lu+Lv;
w=Lv./r;

z=a*w./(a*w+b*(1-w));
t1=betainc(z, a+1, b);
t2=betainc(z ,a, b+1);
A=(1-w).*(1-t1) + w.*t2;

C=exp(-r.*A);