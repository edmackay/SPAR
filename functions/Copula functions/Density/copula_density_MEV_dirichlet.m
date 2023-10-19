function c=copula_density_MEV_dirichlet(u,v,params)

% parameters a,b>0

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
A1=t1+t2-1;

term=gamma(a+b+1)/(gamma(a)*gamma(b));
top=((a*w).^a) .*(b*(1-w)).^b;
bot=(a*w+b*(1-w)).^(a+b+1);
A2=term*top./bot;

% Note that (w*(1-w))^(-1) term in A2 has been cancelled with multiplier

term1=(u.*v).^(A-1);
term2=(A + (1-w).*A1).*(A - w.*A1);
term3=A2./log(u.*v);

c=term1.*(term2 - term3);