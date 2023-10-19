function c=copula_density_Gaussian(u,v,rho)

% parameter 0<a<=1

x=norminv(u);
y=norminv(v);
c = exp(-(rho^2 * (x.^2 + y.^2) - 2*rho*x.*y)/(2*(1-rho^2)))/sqrt(1-rho^2);
c(isnan(c))=0;