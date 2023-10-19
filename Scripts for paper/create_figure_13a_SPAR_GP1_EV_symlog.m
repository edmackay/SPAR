clear
close all

plotsettings_SPAR

a=3; % a>=1
COPULA.name='mev_symlog';
COPULA.params=a; 

MARGINS.X.name='GP'; MARGINS.X.xi=1;
MARGINS.Y.name='GP'; MARGINS.Y.xi=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cartesian space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate joint distribution
rmax=9;
x=logspace(-3,rmax,500);
[X,Y]=meshgrid(x);
fxy = density_xy(X,Y,COPULA,MARGINS);

% calculate SPAR approximation
R=X+Y;
Q=Y./R;
A2=(a - 1).*(Q.*(1-Q)).^(a - 2) .* (Q.^a + (1-Q).^a).^(1/a-2);
fxy_SPAR=R.^(-3).*A2;

% plot comparison
cvals=-40:5:-5;
figure_4panel
contour(log10(x+1),log10(x+1),log10(fxy),cvals)
contour(log10(x+1),log10(x+1),log10(fxy_SPAR),cvals,'k--')
axis([0 1 0 1]*rmax)
xlabel('$x$')
ylabel('$y$')
xticks(0:3:9)
xticklabels({'$10^0$','$10^3$','$10^6$','$10^9$'})
yticks(0:3:9)
yticklabels({'$10^0$','$10^3$','$10^6$','$10^9$'})

hgexport(gcf,['figures/SPAR_GP1_' COPULA.name '_XY_density.eps'], formatEPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Angular-Radial Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate joint distribution
rmax=9;
q=linspace(1e-3,1-1e-3,500);
r=logspace(0,rmax,500);
[Q,R]=meshgrid(q,r);
X=R.*(1-Q);
Y=R.*Q;

frq = density_rq1(R,Q,0,0,COPULA,MARGINS);

% calculate SPAR approximation
A2=(a - 1).*(Q.*(1-Q)).^(a - 2) .* (Q.^a + (1-Q).^a).^(1/a-2);
frq_SPAR=A2.*R.^(-2);

% plot comparison
cvals=-20:2:-2;
figure_4panel
contour(q,log10(r),log10(frq),cvals)
contour(q,log10(r),log10(frq_SPAR),cvals,'k--')
xlabel('$q$')
ylabel('$r$')
yticks(0:3:9)
yticklabels({'$10^0$','$10^3$','$10^6$','$10^9$'})
xlim([0 1])

hgexport(gcf,['figures/SPAR_GP1_' COPULA.name '_RQ_density.eps'], formatEPS);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q0=logspace(-6,-1,50);
q=[0 q0 linspace(0.1,0.9,100) 1-q0 1];
q=sort(q);

figure_4panel
axis([-0.05 1.05 0 3])
xlabel('$q$')
ylabel('$f_Q(q)$')

figure_4panel
axis([-0.05 1.05 0 3])
xlabel('$q$')
ylabel('$\sigma_0(q)$')

ind=0;
for a=[1.4 2 3 4]
    ind=ind+1;
    COPULA.params=a;

    % angular density
    fq=0*q;
    for i=2:length(q)-1
        fq(i)=density_q1(q(i),0,0,COPULA,MARGINS);
    end

    bq=(a-1)*(q.*(1-q)).^(a-2).*(q.^a+(1-q).^a).^(1/a-2);
    sigma=bq./fq;

    figure(3)
    plotcol(q,fq,ind,1,5)

    figure(4)
    plotcol(q,sigma,ind,1,5)
end

figure(3)
L=legend('$\alpha=1.4$', '$\alpha=2$', '$\alpha=3$', '$\alpha=4$');
L.FontSize=12;
L.Box='off';
L.ItemTokenSize = [1 1]*15;
hgexport(gcf,['figures/SPAR_GP1_' COPULA.name '_angular_density.eps'], formatEPS);

figure(4)
hgexport(gcf,['figures/SPAR_GP1_' COPULA.name '_GPscale.eps'], formatEPS);