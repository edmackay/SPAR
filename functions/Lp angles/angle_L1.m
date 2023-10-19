function q1=angle_L1(x,y)

sy=ones(size(x));
sy(y<0)=-1;

r=abs(x)+abs(y);
u=x./r;
q1=sy.*(1-u);
q1(r==0)=0;