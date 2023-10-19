clear

plotsettings_SPAR

x=linspace(-1,1,401);

figure_3panel
grid on; box off
xticks(-1:1)
yticks(-1:1)
axis([-1 1 -1 1]*1.2)

i=0;
colmax=6;
for p=[0.5 1 2 4]
    i=i+1;
    ypos=+(1-abs(x).^p).^(1/p);
    yneg=-(1-abs(x).^p).^(1/p);
    h=plotcol(x,ypos,i,1,colmax);
    h.DisplayName=['$p=' num2str(p) '$'];
    h=plotcol(x,yneg,i,1,7);
    h.HandleVisibility='off';
end
h=plotcol([-1 1 1 -1 -1],[-1 -1 1 1 -1],i+1,1,colmax);
h.DisplayName='$p=\infty$';
xlabel('$x$')
ylabel('$y$')

hgexport(gcf, 'figures/Lp_norm.eps', formatEPS);
