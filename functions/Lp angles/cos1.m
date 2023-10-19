function u=cos1(q)

q=mod(q,4);
q(q>2)=q(q>2)-4;
u=1-abs(q);

