function fxy=density_xy(x,y,COPULA,MARGINS)

fx = margin_pdf(x,MARGINS.X);
fy = margin_pdf(y,MARGINS.Y);
Fx = margin_cdf(x,MARGINS.X);
Fy = margin_cdf(y,MARGINS.Y);
c=copula_density(Fx,Fy,COPULA);
fxy=fx.*fy.*c;
fxy(isnan(fxy))=0;


