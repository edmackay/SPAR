clear

% close all

plotsettings
% set(0,'defaultfigurewindowstyle','docked')
load('FigureFormat_ACP_v3.mat');
formatEPS.Bounds = 'loose';
formatEPS.FontMode='auto';
formatEPS.LineMode='auto';

rho=0.6;
COPULA.name='gaussian';
COPULA.params=rho;

MARGINS.X.name='GP'; MARGINS.X.xi=1;
MARGINS.Y.name='GP'; MARGINS.Y.xi=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cartesian Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate joint distribution
rmax=9;
x=logspace(-3,rmax,500);
[X,Y]=meshgrid(x);
fxy = density_xy(X,Y,COPULA,MARGINS);

% SPAR approximation
R=X+Y;
Q=Y./R;

R(R<1)=NaN;
a=rho/(1+rho);
L = (1-rho^2)^(-0.5) * (4*pi*log(R)).^(-a);
fxy_SPAR = L.* (Q.*(1-Q)).^(a-2) .* R.^(2*a-4);

% plot comparison
cvals=-30:3:-3;
figure_4panel
contour(log10(1+x),log10(1+x),log10(fxy),cvals)
contour(log10(1+x),log10(1+x),log10(fxy_SPAR),cvals,'k--')
axis([0 1 0 1]*rmax)
xlabel('$x$')
ylabel('$y$')
xticks(0:3:9)
xticklabels({'$10^0$','$10^3$','$10^6$','$10^9$'})
yticks(0:3:9)
yticklabels({'$10^0$','$10^3$','$10^6$','$10^9$'})

hgexport(gcf,['figures/SPAR_GP1_' COPULA.name '_XY_density.eps'], formatEPS);

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Angular-Radial Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate joint distribution
rmax=9;
q=linspace(0,1,500);
r=logspace(0,rmax,500);
[Q,R]=meshgrid(q,r);
frq=density_rq1(R,Q,0,0,COPULA,MARGINS);

% calculate SPAR approximation
R(R<1)=NaN;
a=rho/(1+rho);
L = (1-rho^2)^(-0.5) * (4*pi*log(R)).^(-a);
frq_SPAR = L.* (Q.*(1-Q)).^(a-2) .* R.^(2*a-3);

% plot comparison
cvals=-20:2:-2;
figure_4panel
contour(q,log10(r),log10(frq),cvals)
contour(q,log10(r),log10(frq_SPAR),cvals,'k--')
xlabel('$q$')
ylabel('$r$')
yticks(0:3:9)
yticklabels({'$10^0$','$10^3$','$10^6$','$10^9$'})

hgexport(gcf,['figures/SPAR_GP1_' COPULA.name '_RQ_density.eps'], formatEPS);

