function c=copula_density_Nelsen16(u,v,a)

% parameter a>=0
z = u + v - 1 - a*(1./u + 1./v -1);

term1=1 + a./u.^2;
term2=1 + a./v.^2;
term3=(z.^2 + 4*a).^(-3/2);
c=2*a*term1.*term2.*term3; 
