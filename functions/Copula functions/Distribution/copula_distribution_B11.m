function C=copula_distribution_B11(u,v,a)

% parameter 0<=a<=1

C = a*min(u,v)+(1-a)*u*v;