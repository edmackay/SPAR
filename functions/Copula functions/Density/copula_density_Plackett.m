function c=copula_density_Plackett(u,v,a)

% parameter a>=0

b=a-1;
top=a+a*b*(u+v-2*u.*v);
bot=((1+b*(u+v)).^2 - 4*a*b*u.*v).^(3/2);
c=top./bot;
c(isnan(c))=0;
c(isinf(c))=0;