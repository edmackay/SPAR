function Fx = margin_cdf_laplace(x)

Fx = 0.5 + 0.5*sign(x).*(1-exp(-abs(x)));
