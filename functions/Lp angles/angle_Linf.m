function qinf=angle_Linf(x,y)

sy=ones(size(x));
sy(y<0)=-1;

a2 = abs(y)>=abs(x);
a1 = abs(y)<abs(x) & x>0;
a3 = abs(y)<abs(x) & x<0;

r=max(abs(x),abs(y));
u=x./r;
v=y./r;

qinf=zeros(size(x));
qinf(a1)=abs(v(a1))/2;
qinf(a2)=1-u(a2)/2;
qinf(a3)=2-abs(v(a3))/2;

qinf=qinf.*sy;
