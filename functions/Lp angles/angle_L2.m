function q2=angle_L2(x,y)

sy=ones(size(x));
sy(y<0)=-1;

q2=sy.*acos(x./sqrt(x.^2+y.^2))*2/pi;