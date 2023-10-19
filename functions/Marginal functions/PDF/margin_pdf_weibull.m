function fx = margin_pdf_weibull(x,k)

fx=k*x.^(k-1).*exp(-x.^k);
fx(x<0)=0;