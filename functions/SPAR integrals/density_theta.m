function ft=density_theta(theta,x0,y0,COPULA,MARGINS)

% f1=integral(@(r)density_rtheta(r,theta,x0,y0,COPULA,MARGINS),0,15);
% f2=integral(@(r)density_rtheta(r,theta,x0,y0,COPULA,MARGINS),15,inf);
% ft=f1+f2;

ft=integral(@(r)density_rtheta(r,theta,x0,y0,COPULA,MARGINS),0,inf);