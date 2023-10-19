function c=copula_density_student(u,v,params)

rho=params(1);
nu=params(2);

x=tinv(u,nu);
y=tinv(v,nu);
x2=x.^2;
y2=y.^2;

term1=(gamma(nu/2)/gamma((nu+1)/2))^2;
term2=nu/(2*sqrt(1-rho.^2));
term3=(1 + (x2+y2-2*rho.*x.*y)/(nu*(1-rho.^2))).^(-nu/2-1);
term4=((1+x2/nu).*(1+y2/nu)).^((nu+1)/2);
c = term1.*term2.*term3.*term4;
c(isnan(c))=0;