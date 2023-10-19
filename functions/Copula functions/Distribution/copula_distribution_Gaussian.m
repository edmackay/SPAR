function C=copula_distribution_Gaussian(u,v,rho)

% parameter -1<=rho<=1

x=norminv(u);
y=norminv(v);
X=[x(:) y(:)];
R=[1 rho; rho 1];
C = mvncdf(X,[0,0],R);
C = reshape(C,size(u));