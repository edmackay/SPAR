function Frgt=survivor_rcond_L1(rvals,q1,x0,y0,COPULA,MARGINS,fq1)

Frt=0*rvals;
for i=1:length(rvals)
    Frt(i)=integral(@(r)density_rq1(r,q1,x0,y0,COPULA,MARGINS),rvals(i),inf);
end

if nargin<7
    fq1=density_q1(q1,x0,y0,COPULA,MARGINS);
end
Frgt=Frt/fq1;