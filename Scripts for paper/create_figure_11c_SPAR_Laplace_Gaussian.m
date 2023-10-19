clear
close all

plotsettings_SPAR

rho=0.6;
COPULA.name='Gaussian';
COPULA.params=rho; 

MARGINS.X.name='laplace';
MARGINS.Y.name='laplace';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate joint distribution
x=linspace(-20,20,500);
[X,Y]=meshgrid(x);
fxy = density_xy(X,Y,COPULA,MARGINS);

% calculate angular density and threshold
q=linspace(-2,2,401);
fq=0*q;
u=0*q;
zeta=1e-2;
for i=1:length(q)
    disp(num2str(i))
    fq(i)=density_q1(q(i),0,0,COPULA,MARGINS);
    u(i)=quantile_rcond_L1(zeta,q(i),0,0,COPULA,MARGINS,fq(i));
end

% calculate SPAR approximation
R=abs(X)+abs(Y);
Q=angle_L1(X,Y);
W1=cos1(Q);
W2=sin1(Q);

FQ=interp1(q,fq,Q,'pchip');
U=interp1(q,u,Q,'pchip');
R(R<U)=NaN;

LAMBDA=(1 - 2*rho*sign(W1.*W2).*sqrt(abs(W1.*W2)))/(1-rho^2);
SIGMA=1./LAMBDA;
fxy_SPAR=zeta*FQ.*exp(-(R-U).*LAMBDA).*LAMBDA./(SIGMA+U);

% plot comparison
figure_4panel
contour(x,x,log10(fxy),-40:2:-2)
contour(x,x,log10(fxy_SPAR),-40:2:-2,'k--')
set(gca,'YTick',-20:10:20)
set(gca,'xTick',-20:10:20)
title(['$\rho=' num2str(rho) '$, $\zeta=10^{' num2str(log10(zeta)) '}$'] )
xlabel('$x$')
ylabel('$y$')
hgexport(gcf,['figures/SPAR_Laplace_' COPULA.name '_XY_density.eps'], formatEPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate angular density
figure_4panel
ylim([0 0.8])
xlabel('$q$')
ylabel('$f_Q(q)$')

q=linspace(-2,2,401);
ind=0;
for rho=0:0.2:0.8
    ind=ind+1;
    COPULA.params=rho; 
    fq=0*q;
    for i=1:length(q)
        disp(num2str(i))
        fq(i)=density_q1(q(i),0,0,COPULA,MARGINS);
    end
    h=plotcol(q,fq,ind,1,6);
    h.DisplayName=['$\rho=' num2str(rho) '$'];
end
L=legend;
L.Box='off';
L.Location="northeast";
L.ItemTokenSize = [1 1]*15;
L.Position=[0.3233    0.5792    0.2547    0.3380];
hgexport(gcf,['figures/SPAR_Laplace_' COPULA.name '_angular_density.eps'], formatEPS);

% calculate limit sets and rate parameter
figure_4panel
axis([-2 2 0 2])
xlabel('$q$')
ylabel('$\sigma(q)$')

figure_4panel
axis([-1 1 -1 1])
xlabel('$x$')
ylabel('$y$')
xticks(-1:1)
yticks(-1:1)

q=linspace(-2,2,401);
w1=cos1(q);
w2=sin1(q);

ind=0;
for rho=0:0.2:0.8
    ind=ind+1;

    lambda=(1 - 2*rho*sign(w1.*w2).*sqrt(abs(w1.*w2)))/(1-rho^2);

    figure(3)
    plotcol(q,1./lambda,ind,1,6);
    figure(4)
    plotcol(cos1(q)./lambda,sin1(q)./lambda,ind,1,6);
end

figure(3)
hgexport(gcf,['figures/SPAR_Laplace_' COPULA.name '_GPscale.eps'], formatEPS);
figure(4)
hgexport(gcf,['figures/SPAR_Laplace_' COPULA.name '_limitset.eps'], formatEPS);