function fq=density_q1(q1,x0,y0,COPULA,MARGINS)

% f1=integral(@(r)density_rq1(r,q1,x0,y0,COPULA,MARGINS),0,15);
% f2=integral(@(r)density_rq1(r,q1,x0,y0,COPULA,MARGINS),15,inf);
% ft=f1+f2;

fq=integral(@(r)density_rq1(r,q1,x0,y0,COPULA,MARGINS),0,inf);