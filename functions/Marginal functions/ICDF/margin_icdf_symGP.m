function x = margin_icdf_symGP(p,xi)

x = sign(p-0.5).*gpinv(2*abs(p-0.5),xi,1,0);
