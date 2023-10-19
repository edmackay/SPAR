function Frgt=exceedance_rcond(rvals,theta,x0,y0,rmax,COPULA,MARGINS,ftheta)

Frt=0*rvals;
for i=1:length(rvals)
    Frt(i)=integral(@(r)density_rtheta(r,theta,x0,y0,COPULA,MARGINS),rvals(i),rmax);
end

if nargin<8
    ftheta=density_theta(theta,x0,y0,COPULA,MARGINS);
end
Frgt=Frt/ftheta;