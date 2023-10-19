clear
close all
plotsettings_SPAR

%%%%%%%%%%%%%%%%%%%%%%%%%% Comparison between SPAR and true density %%%%%%%%%%%%%%%%%%%%%%%%%%
% distribution parameters
nu=2;
rho=0.6;
C=[1 rho; rho 1];

% calculate joint distribution
x=linspace(-30,30,500);
[X,Y]=meshgrid(x);
fxy=mvtpdf([X(:) Y(:)],C,nu);
fxy=reshape(fxy,size(X));

% calculate SPAR approximation
R=sqrt(X.^2+Y.^2);
THETA=atan2(Y,X);

% angular density
FT=sqrt(1-rho^2)./(2*pi*(1-rho*sin(2*THETA))); 

% GP parameters
zeta=0.05;
ALPHA=(1-rho^2)./(1-rho*sin(2*THETA));
U=sqrt(nu*ALPHA*(zeta^(-2/nu)-1));
SIGMA=U/nu;
XI=1/nu;
R(R<U)=NaN;

% SPAR approximation
GP=gppdf(R,XI,SIGMA,U);
fxy_SPAR=zeta.*FT.*GP./R;

% plot comparison
figure_3panel
contour(x,x,log10(fxy),-7:0)
contour(x,x,log10(fxy_SPAR),-7:0,'k--')
set(gca,'YTick',-30:15:30)
set(gca,'xTick',-30:15:30)
caxis([-7 0])
xlabel('$x$')
ylabel('$y$')
hgexport(gcf, 'figures/SPAR_student_fxy.eps', formatEPS);

%%%%%%%%%%%%%%%%%%%%%%%%%% Angular density and GP scale parameter %%%%%%%%%%%%%%%%%%%%%%%%%%
theta=linspace(-pi,pi,1000);
ft=sqrt(1-rho^2)./(2*pi*(1-rho*sin(2*theta))); 
alpha=(1-rho^2)./(1-rho*sin(2*theta));
u=sqrt(nu*alpha*(zeta^(-2/nu)-1));
sigma=u/nu;

figure_3panel
plot(theta,ft)
xlim([-pi pi])
xticks(-pi:pi/2:pi)
xticklabels({'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'})
axis square
xlabel('$\theta$')
ylabel('$f_{\Theta}(\theta)$')
set(gca,'XTickLabelRotation',0)
ylim([0 0.4])
hgexport(gcf, 'figures/SPAR_student_ft.eps', formatEPS);

figure_3panel
plot(theta,sigma)
xlim([-pi pi])
xticks(-pi:pi/2:pi)
xticklabels({'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'})
xlabel('$\theta$')
ylabel('$\sigma(\theta)$')
set(gca,'XTickLabelRotation',0)
ylim([0 5])
hgexport(gcf,'figures/SPAR_student_sigma.eps', formatEPS);

