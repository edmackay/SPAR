function Fx = margin_cdf_symGP(x,xi)

Fx = 0.5 + 0.5*sign(x).*gpcdf(abs(x),xi,1,0);
