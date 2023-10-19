function c=copula_density_Clayton(u,v,a)

% parameter a>0

z=u.^(-a)+v.^(-a)-1;
c=(a+1) * ((u.*v).^(-a-1)).* (z.^(-2-1/a));
c(z<0)=0;
c(isnan(c))=0;
c(isinf(c))=0;