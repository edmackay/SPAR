function frq=density_rq1(r,q1,x0,y0,COPULA,MARGINS)

x=x0+r.*cos1(q1);
y=y0+r.*sin1(q1);

fxy=density_xy(x,y,COPULA,MARGINS);
frq=r.*fxy;

