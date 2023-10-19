clear

plotsettings_SPAR

a=3;
COPULA.name='Nelsen15';
COPULA.params=a;
COPULA.lowertail=1;

MARGINS.X.name='GP'; MARGINS.X.xi=-1/a;
MARGINS.Y.name='GP'; MARGINS.Y.xi=-1/a;

%%%%%%%%%%%%%%%%%% Copula contour plot %%%%%%%%%%%%%%%%%% 
u=linspace(0,1,1000);
[U,V]=meshgrid(u);
C=copula_distribution_Nelsen15(U,V,a);

figure_3panel
contour(u,u,C,0.1:0.1:1)
contour(u,u,C,[1 1]*1e-3,'k','LineWidth',2)
xlabel('$u$')
ylabel('$v$')

hgexport(gcf, 'figures/SPAR_GPneg_Nelsen15_copula.eps', formatEPS);

%%%%%%%%%%%%%%%%%% Cartesian density plot %%%%%%%%%%%%%%%%%% 
x=linspace(0,a,500);
[X,Y]=meshgrid(x);

fxy=density_xy(X,Y,COPULA,MARGINS);

q=linspace(0,1,1000);
rF=a.*(q.^a+(1-q).^a).^(-1/a);

figure_3panel
contour(x,x,log10(fxy),-6:0.5:0)
plot(rF.*(1-q),rF.*q,'k-','LineWidth',2)
xlabel('$x$')
ylabel('$y$')

hgexport(gcf, 'figures/SPAR_GPneg_Nelsen15_XYdens.eps', formatEPS);

%%%%%%%%%%%%%%%%%%%%%%%% angular density %%%%%%%%%%%%%%%%%%%%%%%% 
q=linspace(0,1,1000);

figure_3panel
xlabel('$q$')
ylabel('$f_Q(q)$')
axis([0 1 0 5])
for a=1:4
    fq=a*(q.*(1-q)).^(a-1) .* ((1-q).^a + q.^a).^(-2);

    plotcol(q,fq,a,1,5)
end
L=legend('$\alpha=1$','$\alpha=2$','$\alpha=3$','$\alpha=4$');
L.Box='off';
L.Location='northwest';
% L.FontSize=12;

hgexport(gcf, 'figures/SPAR_GPneg_Nelsen15_fQ.eps', formatEPS);
