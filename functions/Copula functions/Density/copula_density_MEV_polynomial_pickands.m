function c=copula_density_MEV_polynomial_pickands(u,v,params)

a=params(1);
b=params(2);

if (a<0) || (a+3*b<0) || (a+b>1) || (a+2*b>1)
    error('parameters out of range')
end

w=log(v)./log(u.*v);

A = 1 - (a+b)*w + a*w.^2 + b*w.^3;
A1 = -(a+b) + 2*a*w + 3*b*w.^2;
A2 = 2*a + 6*b*w;

term1=(u.*v).^(A-1);

term2=(A + (1-w).*A1).*(A - w.*A1);

term3=w.*(1-w)./log(u.*v);

c=term1.*(term2 - term3.*A2);

c(isnan(c))=0;
c(isinf(c))=0;

