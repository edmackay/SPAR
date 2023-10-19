function Frgt=distribution_rcond_L1(rvals,q1,x0,y0,COPULA,MARGINS,fq1)

Frt=0*rvals;
for i=1:length(rvals)
    if rvals(i)>15
        f1=integral(@(r)density_rq1(r,q1,x0,y0,COPULA,MARGINS),0,15);
        f2=integral(@(r)density_rq1(r,q1,x0,y0,COPULA,MARGINS),15,rvals(i));
        Frt(i)=f1+f2;
    else
        Frt(i)=integral(@(r)density_rq1(r,q1,x0,y0,COPULA,MARGINS),0,rvals(i));
    end
end

if nargin<7
    fq1=density_q1(q1,x0,y0,COPULA,MARGINS);
end
Frgt=Frt/fq1;