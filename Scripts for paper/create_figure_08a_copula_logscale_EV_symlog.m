clear
close all

addpath(genpath('..\functions'))
plotsettings_SPAR

a=2;
u=logspace(-10,0,200);
[U,V]=meshgrid(u);
T=U+V;
W=V./T;
R=U.*V;
Q=log(V)./log(R);

%%%%%%%%%%%%%%%%%%%%%% lower tail %%%%%%%%%%%%%%%%%%%%%%
Lu=-log(U);
Lv=-log(V);
z=Lu.^a + Lv.^a;
C=exp(-z.^(1/a));
kappa=2^(1/a);
C_approx_HJ=T.^kappa .* (W.*(1-W)).^(kappa/2);

figure_3panel
contour(log10(u),log10(u),log10(C),-15:0)
contour(log10(u),log10(u),log10(C_approx_HJ),-15:0,'k--')
xlabel('$u$')
ylabel('$v$')
xticks(-10:5:0)
xticklabels({'$10^{-10}$','$10^{-5}$','$10^{0}$'})
yticks(-10:5:0)
yticklabels({'$10^{-10}$','$10^{-5}$','$10^{0}$'})

hgexport(gcf, 'figures/Copula_EV_symlog_lowertail_HJ.eps', formatEPS);

%%%%%%%%%%%%%%%%%%%%%% upper tail %%%%%%%%%%%%%%%%%%%%%%
Lu=-log(1-U);
Lv=-log(1-V);
z=Lu.^a + Lv.^a;
C11 = T - 1 + exp(-z.^(1/a));
C11_approx_WT=R.^(max(Q,1-Q));

figure_3panel
contour(log10(u),log10(u),log10(C11),-15:0)
contour(log10(u),log10(u),log10(C11_approx_WT),-15:0,'k--')
xlabel('$u$')
ylabel('$v$')
xticks(-10:5:0)
xticklabels({'$10^{-10}$','$10^{-5}$','$10^{0}$'})
yticks(-10:5:0)
yticklabels({'$10^{-10}$','$10^{-5}$','$10^{0}$'})

hgexport(gcf, 'figures/Copula_EV_symlog_uppertail_WT.eps', formatEPS);