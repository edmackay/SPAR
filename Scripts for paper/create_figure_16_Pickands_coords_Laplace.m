clear

plotsettings_SPAR

FS=14;

figure('Position',[100 100 1200 300])
subplot(1,2,1)
hold on; grid on; zoom on
set(gca,'yscale','log')
set(gca,'yminortick','off')
set(gca,'yminorgrid','off')
ylim([1e-6 1e2])
yticks(logspace(-8,2,6))
xlabel('$q$')
ylabel('$\tau$')
set(gca,'fontsize',FS)

subplot(1,2,2)
hold on; grid on; zoom on
xlabel('$q$')
ylabel('$\psi$')
set(gca,'fontsize',FS)

for r=5:5:25
    q=-2:1e-3:2;
    u=margin_cdf_laplace(r*cos1(q));
    v=margin_cdf_laplace(r*sin1(q));
    t=-log(u.*v);
    w=-log(v)./t;
    subplot(1,2,1)
    plotcol(q,t,r,5,30)
    subplot(1,2,2)
    h=plotcol(q,w,r,5,30);
    h.DisplayName=['$r=' num2str(r) '$'];
end
L=legend;
L.Box='off';
L.FontSize=FS;

hgexport(gcf, 'figures/Pickands_coords.eps', formatEPS);

