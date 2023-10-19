clear

plotsettings_SPAR

a=-0.2;
COPULA.name='Clayton';
COPULA.params=a;
COPULA.lowertail=1;

MARGINS.X.name='GP'; MARGINS.X.xi=a;
MARGINS.Y.name='GP'; MARGINS.Y.xi=a;

% %%%%%%%%%%%%%%%%%% Copula contour plot %%%%%%%%%%%%%%%%%% 
% u=linspace(0,1,1000);
% [U,V]=meshgrid(u);
% C=copula_distribution_Clayton(U,V,a);
% 
% figure_3panel
% contour(u,u,C,0.1:0.1:1)
% contour(u,u,C,[1 1]*1e-3,'k','LineWidth',2)
% xlabel('$u$')
% ylabel('$v$')
% 
% hgexport(gcf, 'figures/SPAR_GPneg_Clayton_copula.eps', formatEPS);
% 
% %%%%%%%%%%%%%%%%%% Cartesian density plot %%%%%%%%%%%%%%%%%% 
% x=linspace(0,-1/a,500);
% [X,Y]=meshgrid(x);
% 
% fxy=density_xy(X,Y,COPULA,MARGINS);
% 
% figure_3panel
% contour(x,x,log10(fxy),-6:0.5:0)
% plot([-1/a 0],[0 -1/a],'k-','LineWidth',2)
% xlabel('$x$')
% ylabel('$y$')
% 
% hgexport(gcf, 'figures/SPAR_GPneg_Clayton_XYdens.eps', formatEPS);
% 
% %%%%%%%%%%%%%%%%%% Density vs GP %%%%%%%%%%%%%%%%%% 
% 
% rF=-1/a;
% s=logspace(-2,log10(rF),1000);
% r=rF-s;
% 
% q=0.5;
% frq=density_rq1(r,q,0,0,COPULA,MARGINS);
% 
% xi=a/(1+a);
% sigma=(-xi)^(-xi)*(1-xi)^(xi+1);
% fgp=sigma^(1/xi) * (-xi*s).^(-1-1/xi);
% 
% figure_3panel
% set(gca,'yScale','log')
% set(gca,'xScale','log')
% plot(s/rF,frq)
% plot(s/rF,fgp,'k--')
% axis([1e-2 1 1e-6 1])
% xlabel('$1-r/r_F$')
% L=legend('$f_{R|Q}(r|q)$','$f_{GP}(r; \xi, \sigma)$');
% L.Box='off';
% L.Location='northwest';
% L.FontSize=14;
% hgexport(gcf, 'figures/SPAR_GPneg_Clayton_RQdens.eps', formatEPS);
% 
% %%%%%%%%%%%%%%%%%%%%%%%% angular density %%%%%%%%%%%%%%%%%%%%%%%% 

figure_3panel
xlabel('$q$')
ylabel('$f_Q(q)$')
axis([0 1 0 2])
plot([0 1],[1 1],'k')
L=legend('$\alpha\in[-1,0)$');
L.Box='off';
L.Location='northwest';
% L.FontSize=12;

hgexport(gcf, 'figures/SPAR_GPneg_Clayton_fQ.eps', formatEPS);