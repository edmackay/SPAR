clear

plotsettings_SPAR

a=10; % -inf<=a<=inf
COPULA.name='Frank';
COPULA.params=a; 

MARGINS.X.name='laplace';
MARGINS.Y.name='laplace';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% calculate joint distribution
x=linspace(-20,20,500);
[X,Y]=meshgrid(x);
fxy = density_xy(X,Y,COPULA,MARGINS);

% calculate angular density and threshold
q=linspace(-2,2,401);
fq=0*q;
u=0*q;
zeta=1e-4;
disp('Calculating angular density and threshold for angle...')
for i=1:length(q)
    disp(num2str(i))
    fq(i)=density_q1(q(i),0,0,COPULA,MARGINS);
    u(i)=quantile_rcond_L1(zeta,q(i),0,0,COPULA,MARGINS,fq(i));
end

% calculate SPAR approximation
R=abs(X)+abs(Y);
Q=angle_L1(X,Y);
FQ=interp1(q,fq,Q,'pchip');
U=interp1(q,u,Q,'pchip');
R(R<U)=NaN;
fxy_SPAR=zeta*FQ.*exp(U-R)./(1+U);

% plot comparison
figure_4panel
contour(x,x,log10(fxy),-20:2:-2)
contour(x,x,log10(fxy_SPAR),-20:2:-2,'k--')
set(gca,'YTick',-20:10:20)
set(gca,'xTick',-20:10:20)
title(['$\alpha=' num2str(a) '$, $\zeta=10^{' num2str(log10(zeta)) '}$'] )
xlabel('$x$')
ylabel('$y$')
hgexport(gcf,['figures/SPAR_Laplace_' COPULA.name '_XY_density.eps'], formatEPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate angular density
figure_4panel
ylim([0 1])
xlabel('$q$')
ylabel('$f_Q(q)$')

q=linspace(-2,2,401);
disp('Calculating angular density for angle...')
for a=10:-2:2
    COPULA.params=a; 
    fq=0*q;
    for i=1:length(q)
        disp(num2str(i))
        fq(i)=density_q1(q(i),0,0,COPULA,MARGINS);
    end
    h=plotcol(q,fq,a,2,12);
    h.DisplayName=['$\alpha=' num2str(a) '$'];
end
L=legend;
L.Box='off';
L.Position=[0.3136    0.5837    0.3080    0.3380];
hgexport(gcf,['figures/SPAR_Laplace_' COPULA.name '_angular_density.eps'], formatEPS);

% calculate limit sets and rate parameter
q=linspace(-2,2,401);
lambda=1+0*q;

figure_4panel
axis([-2 2 0 2])
xlabel('$q$')
ylabel('$\sigma(q)$')
plotcol(q,lambda,1,1,5);
hgexport(gcf,['figures/SPAR_Laplace_' COPULA.name '_GPscale.eps'], formatEPS);

figure_4panel
axis([-1 1 -1 1])
xticks(-1:1)
yticks(-1:1)
xlabel('$x$')
ylabel('$y$')
plotcol(cos1(q)./lambda,sin1(q)./lambda,1,1,5);
hgexport(gcf,['figures/SPAR_Laplace_' COPULA.name '_limitset.eps'], formatEPS);