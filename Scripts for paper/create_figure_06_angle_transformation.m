clear
close all
addpath('..')

plotsettings_SPAR

figure_3panel
plot([1 0 -1 0 1],[0 1 0 -1 0],'k--','LineWidth',1)
plot([0.5 0.5 0 -0.5 -0.5 0 0.5],[-0.5 0.5 1 0.5 -0.5 -1 -0.5],'k','LineWidth',2)
axis([-1 1 -1 1]*1.2)
axis square
grid on
xlabel('$x$')
ylabel('$y$')
L=legend('$\mathcal{U}_1$','$\mathcal{S}_a$');
L.Box='off';
L.FontSize=14;
L.Position=[0.6891    0.7287    0.1594    0.1174];
% hgexport(gcf, '/figures/badcurve.eps', formatEPS);

a=(1+sqrt(2));
q=linspace(0,1/a,1000);
dq1dq=a.*(1+a*q).^-2;

figure_3panel
hold on
plot([0 1],[1 1],'k--')
plot(q,dq1dq,'k','LineWidth',2)
plot([1/a 1],[1 1]*a/(2*sqrt(2)),'k','LineWidth',2,'HandleVisibility','off')
xlabel('$q$')
ylim([0 2.5])
L=legend('$f_{Q_1}(q)$','$f_{Q_a}(q)$');
L.Box='off';
L.FontSize=14;

% hgexport(gcf, '/figures/badcurve_density.eps', formatEPS);
