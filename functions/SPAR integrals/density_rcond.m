function frgt=density_rcond(r,theta,x0,y0,COPULA,MARGINS,ft)

if nargin<7
    ft=density_theta(theta,x0,y0,COPULA,MARGINS);
end
frt=density_rtheta(r,theta,x0,y0,COPULA,MARGINS);
frgt=frt/ft;