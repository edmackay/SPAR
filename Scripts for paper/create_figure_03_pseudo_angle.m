clear
close all

plotsettings_SPAR

x=linspace(-1,1,1000);
p=3;
ypos=+(1-abs(x).^p).^(1/p);
yneg=-(1-abs(x).^p).^(1/p);

amin=-1.2;
amax=1.2;

x0=0.75;
y0=(1-abs(x0).^p).^(1/p);

FontSize=16;

figure('Position',[100 100 500 450])
hold on
axis square
axis([amin amax amin amax])
axis off
plot(x(x>x0),ypos(x>x0),'k','LineWidth',2)
plot(x(x<=x0),ypos(x<=x0),'k--')
plot(x,yneg,'k--')
plot([0 0],[amin amax],'k')
plot([amin amax],[0 0],'k')
plot([x0 x0],[0 y0],'r-','LineWidth',2)
plot([0 x0],[0 0],'b-','LineWidth',2)
plot([0 x0],[0 y0],'k-','LineWidth',0.5)
scatter(x0,y0,'ko','filled')
annotation('arrow',[0.228 0.288],[0.89 0.806]);
annotation('arrow',[0.771 0.750],[0.770 0.801],'HeadStyle','deltoid');
text(-1.475,1.18,'Curve $\mathcal{S}_*$','FontSize',FontSize)
text(0.741,0.936,'$\mathbf{w}=\mathbf{x}/\mathcal{R}_*(\mathbf{x})$','FontSize',FontSize)
text(1.04, 0.459,'$\ell_*(\mathbf{w})$','FontSize',FontSize)
text(0.185,-0.107,'$\cos_*(q)$','FontSize',FontSize,'Color','b')
text(0.64,0.147,'$\sin_*(q)$','FontSize',FontSize,'Color','r','Rotation',90)
text(-0.829375,-0.379,'$$q=\frac{4}{\mathcal{C}_*}\ell_*(\mathbf{w})$$','FontSize',FontSize)

hgexport(gcf, 'figures/pseudo_angle.eps', formatEPS);
