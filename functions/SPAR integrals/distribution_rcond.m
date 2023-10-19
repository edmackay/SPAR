function Frgt=distribution_rcond(rvals,theta,x0,y0,COPULA,MARGINS,ftheta)

Frt=0*rvals;
for i=1:length(rvals)
    if rvals(i)>15
        f1=integral(@(r)density_rtheta(r,theta,x0,y0,COPULA,MARGINS),0,15);
        f2=integral(@(r)density_rtheta(r,theta,x0,y0,COPULA,MARGINS),15,rvals(i));
        Frt(i)=f1+f2;
    else
        Frt(i)=integral(@(r)density_rtheta(r,theta,x0,y0,COPULA,MARGINS),0,rvals(i));
    end
end

if nargin<7
    ftheta=density_theta(theta,x0,y0,COPULA,MARGINS);
end
Frgt=Frt/ftheta;