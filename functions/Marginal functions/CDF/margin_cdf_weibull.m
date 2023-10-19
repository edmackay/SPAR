function Fx = margin_cdf_weibull(x,k)

Fx = 1-exp(-x.^k);
Fx(x<0)=0;
