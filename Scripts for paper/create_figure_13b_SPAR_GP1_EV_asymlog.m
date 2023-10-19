clear
close all

plotsettings_SPAR

a=5;
psi1=0.9;
psi2=0.1;
COPULA.name='MEV_asymlog';
COPULA.params=[a,psi1,psi2];

MARGINS.X.name='GP'; MARGINS.X.xi=1;
MARGINS.Y.name='GP'; MARGINS.Y.xi=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cartesian space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate joint distribution
rmax=15;
x=logspace(-2,rmax,500);
[X,Y]=meshgrid(x);
fxy = density_xy(X,Y,COPULA,MARGINS);

% calculate SPAR approximation
R=abs(X)+abs(Y);
W=X./R;
Z = (psi1*(1-W)).^a + (psi2*W).^a;
A2 = (a-1) * (psi1*psi2)^a * (W.*(1-W)).^(a-2) .* Z.^(1/a-2);
fxy_SPAR=A2.*R.^(-3);

% plot comparison
cvals=-50:5:-5;
figure_4panel
contour(log10(1+x),log10(1+x),log10(fxy),cvals)
contour(log10(1+x),log10(1+x),log10(fxy_SPAR),cvals,'k--')
axis([0 1 0 1]*rmax)
xlabel('$x$')
ylabel('$y$')
xticklabels({'$10^0$','$10^5$','$10^{10}$','$10^{15}$'})
yticks(0:5:15)
yticklabels({'$10^0$','$10^5$','$10^{10}$','$10^{15}$'})

hgexport(gcf,['figures/SPAR_GP1_' COPULA.name '_XY_density.eps'], formatEPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Angular-Radial Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate joint distribution
rmax=15;
q=linspace(1e-4,1-1e-4,500);
r=logspace(0,rmax,500);
[Q,R]=meshgrid(q,r);

frq = density_rq1(R,Q,0,0,COPULA,MARGINS);

% calculate SPAR approximation
W=1-Q;
Z = (psi1*(1-W)).^a + (psi2*W).^a;
A2 = (a-1) * (psi1*psi2)^a * (W.*(1-W)).^(a-2) .* Z.^(1/a-2);
frq_SPAR=A2.*R.^(-2);

% plot comparison
cvals=-35:5:-5;
figure_4panel
contour(q,log10(r),log10(frq),cvals)
contour(q,log10(r),log10(frq_SPAR),cvals,'k--')
xlabel('$q$')
ylabel('$r$')
yticks(0:5:15)
yticklabels({'$10^0$','$10^5$','$10^{10}$','$10^{15}$'})
xlim([0 1])

hgexport(gcf,['figures/SPAR_GP1_' COPULA.name '_RQ_density.eps'], formatEPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q0=logspace(-10,-1,400);
q=[0 q0 linspace(0.1,0.9,200) 1-q0 1];
q=sort(q);

figure_4panel
axis([-0.05 1.05 0 3])
xlabel('$q$')
ylabel('$f_Q(q)$')

figure_4panel
axis([-0.05 1.05 0 8])
xlabel('$q$')
ylabel('$\sigma_0(q)$')

ind=0;
for a=[1.5 2 3 4 5]
    ind=ind+1;
    COPULA.params=[a psi1 psi2];

    % angular density
    fq=0*q;
    for i=2:length(q)-1
        fq(i)=density_q1(q(i),0,0,COPULA,MARGINS);
    end

    w=1-q;
    z = (psi1*(1-w)).^a + (psi2*w).^a;
    bq = (a-1) * (psi1*psi2)^a * (w.*(1-w)).^(a-2) .* z.^(1/a-2);
    sigma=bq./fq;
    sigma(1)=0;
    sigma(end)=0;

    figure(3)
    plotcol(q,fq,ind,1,6)

    figure(4)
    plotcol(q,sigma,ind,1,6)
end

figure(3)
L=legend('$\alpha=1.5$', '$\alpha=2$', '$\alpha=3$', '$\alpha=4$', '$\alpha=5$');
L.FontSize=12;
L.Box='off';
L.Position=[0.5086    0.5863    0.3142    0.3228];
hgexport(gcf,['figures/SPAR_GP1_' COPULA.name '_angular_density.eps'], formatEPS);

figure(4)
hgexport(gcf,['figures/SPAR_GP1_' COPULA.name '_GPscale.eps'], formatEPS);

