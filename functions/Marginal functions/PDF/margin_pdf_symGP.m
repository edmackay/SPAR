function fx = margin_pdf_symGP(x,xi)

fx = 0.5*gppdf(abs(x),xi,1,0);
