function C=copula_distribution_Raftery(u,v,a)

% parameter 0<=a<=1

x=min(u,v);
y=max(u,v);
C = x - ((1-a)/(1+a))*(x.^(1/(1-a))).*(y.^(-a/(1-a)) - y.^(1/1-a));