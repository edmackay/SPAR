function x = margin_icdf_laplace(p)

x = -sign(p-0.5).*log(1-2*abs(p-0.5));