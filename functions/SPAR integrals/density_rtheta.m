function frt=density_rtheta(r,theta,x0,y0,COPULA,MARGINS)

x=x0+r.*cos(theta);
y=y0+r.*sin(theta);

fxy=density_xy(x,y,COPULA,MARGINS);
frt=r.*fxy;

