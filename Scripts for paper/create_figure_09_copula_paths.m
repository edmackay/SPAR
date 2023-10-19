clear
close all

addpath(genpath('..\functions'))
plotsettings_SPAR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exponential 
figure_3panel
axis([0 1 0 1])
xlabel('$u$')
ylabel('$v$')

r=linspace(0,50,1000);
for q=0.1:0.1:0.9
    u=1-exp(-r*(1-q));
    v=1-exp(-r*q);
    plot(u,v,'k')
end

hgexport(gcf, 'figures/Copula_paths_exponential.eps', formatEPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Laplace 
figure_3panel
axis([0 1 0 1])
xlabel('$u$')
ylabel('$v$')

r=linspace(0,50,1000);
for q=-2:0.1:2
    x=r*cos1(q);
    y=r*sin1(q);
    u=margin_cdf_laplace(x);
    v=margin_cdf_laplace(y);
    plot(u,v,'k')
end

hgexport(gcf, 'figures/Copula_paths_Laplace.eps', formatEPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP1
figure_3panel
axis([0 1 0 1])
xlabel('$u$')
ylabel('$v$')

r=logspace(-3,4,1000);
for q=0.1:0.1:0.9
    u=1-(1+r*(1-q)).^-1;
    v=1-(1+r*q).^-1;
    plot(u,v,'k')
end

hgexport(gcf, 'figures/Copula_paths_GP1.eps', formatEPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pareto (0,0)
figure_3panel
axis([0 1 0 1])
xlabel('$u$')
ylabel('$v$')

for q=0.1:0.1:0.9
    rmin=max(1/q,1/(1-q));
    r=logspace(log10(rmin),4,1000);
    u=1-(r*(1-q)).^-1;
    v=1-(r*q).^-1;
    plot(u,v,'k')
end

hgexport(gcf, 'figures/Copula_paths_pareto00.eps', formatEPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP negative
figure_3panel
axis([0 1 0 1])
xlabel('$u$')
ylabel('$v$')

xi=-1;

for q=0.1:0.1:0.9
    rmax=min(-1/(xi*(1-q)),-1/(xi*q));
    r=linspace(0,rmax,1000);
    u=gpcdf(r*(1-q),xi,1,0);
    v=gpcdf(r*q,xi,1,0);
    plot(u,v,'k')
end

hgexport(gcf, 'figures/Copula_paths_GPm1.eps', formatEPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SGP
figure_3panel
axis([0 1 0 1])
xlabel('$u$')
ylabel('$v$')

xi=1;

for q=-2:0.1:2
    r=logspace(-3,10,1000);
    x=r*cos1(q);
    y=r*sin1(q);
    u=margin_cdf_symGP(x,xi);
    v=margin_cdf_symGP(y,xi);
    plot(u,v,'k')
end

hgexport(gcf, 'figures/Copula_paths_SGP1.eps', formatEPS);