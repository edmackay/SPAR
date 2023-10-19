clear

plotsettings_SPAR

COPULA.name='Gaussian';
COPULA.params=0.6; % rho in [-1,1]

% rho=0.2;
% nu=2;
% COPULA.name='student';
% COPULA.params=[rho nu];

% COPULA.name='MEV_symlog';
% COPULA.params=2; % a>=1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=linspace(0,1,500);
[U,V]=meshgrid(u);
c = copula_density(U,V,COPULA);
cs=c./(1+c);

figure_3panel
contour(u,u,cs,0:0.1:1)
axis([0 1 0 1])
caxis([0 1.1])
xlabel('$u$')
ylabel('$v$')

hgexport(gcf,['figures/copula_density_' COPULA.name '.eps'], formatEPS);