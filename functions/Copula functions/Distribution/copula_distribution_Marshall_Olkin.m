function C=copula_distribution_Marshall_Olkin(u,v,a,b)

% parameter a must be in [0,1]

C=min(u.^a,v.^b).*(u.^(1-a)).*(v.^(1-b));