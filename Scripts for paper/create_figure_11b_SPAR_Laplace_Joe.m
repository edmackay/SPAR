clear
close all

plotsettings_SPAR

a=3; % a>0
COPULA.name='Joe';
COPULA.params=a; 

MARGINS.X.name='laplace';
MARGINS.Y.name='laplace';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% calculate joint distribution
x=linspace(-20,20,500);
[X,Y]=meshgrid(x);
fxy = density_xy(X,Y,COPULA,MARGINS);

% calculate angular density and threshold
q=linspace(-2,2,401);
fq=0*q;
u=0*q;
zeta=0.05;
for i=1:length(q)
    disp(num2str(i))
    fq(i)=density_q1(q(i),0,0,COPULA,MARGINS);
    u(i)=quantile_rcond_L1(zeta,q(i),0,0,COPULA,MARGINS,fq(i));
end

% calculate SPAR approximation
R=abs(X)+abs(Y);
Q=angle_L1(X,Y);
C=abs(cos1(Q));
S=abs(sin1(Q));

FQ=interp1(q,fq,Q,'pchip');
U=interp1(q,u,Q,'pchip');
R(R<U)=NaN;

Q3 = Q>-2 & Q<=-1;
Q4 = Q>-1 & Q<= 0;
Q1 = Q> 0 & Q<= 1;
Q2 = Q> 1 & Q<= 2;

LAMBDA=0*C;
LAMBDA(Q1)=a+(1-2*a)*min(C(Q1),S(Q1));
LAMBDA(Q2)=1+(a-1)*S(Q2);
LAMBDA(Q3)=1;
LAMBDA(Q4)=1+(a-1)*C(Q4);
SIGMA=1./LAMBDA;
fxy_SPAR=zeta*FQ.*exp(-(R-U).*LAMBDA).*LAMBDA./(SIGMA+U);

% plot comparison
figure_4panel
contour(x,x,log10(fxy),-34:2:-2)
contour(x,x,log10(fxy_SPAR),-34:2:-2,'k--')
set(gca,'YTick',-20:10:20)
set(gca,'xTick',-20:10:20)
title(['$\alpha=' num2str(a) '$, $\zeta=' num2str(zeta) '$'] )
xlabel('$x$')
ylabel('$y$')
hgexport(gcf,['figures/SPAR_Laplace_' COPULA.name '_XY_density.eps'], formatEPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate angular density
figure_4panel
ylim([0 1.5])
xlabel('$q$')
ylabel('$f_Q(q)$')

q=linspace(-2,2,401);
ind=0;
for a=[1 1.5 2 3]
    ind=ind+1;
    COPULA.params=a; 
    fq=0*q;
    for i=1:length(q)
        disp(num2str(i))
        fq(i)=density_q1(q(i),0,0,COPULA,MARGINS);
    end
    h=plotcol(q,fq,ind,1,5);
    h.DisplayName=['$\alpha=' num2str(a) '$'];
end
L=legend;
L.Box='off';
L.Location="northwest";
hgexport(gcf,['figures/SPAR_Laplace_' COPULA.name '_angular_density.eps'], formatEPS);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
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
q3 = q>=-2 & q<=-1;
q4 = q>-1 & q<= 0;
q1 = q> 0 & q<= 1;
q2 = q> 1 & q<= 2;
c=abs(cos1(q));
s=abs(sin1(q));

ind=0;
for a=[1 1.5 2 3]
    ind=ind+1;

    lambda=0*q;
    if a==1
        lambda=1+0*q;
    else
        lambda(q1)=a+(1-2*a)*min(c(q1),s(q1));
        lambda(q2)=c(q2)+a*s(q2);
        lambda(q3)=1;
        lambda(q4)=a*c(q4)+s(q4);
    end

    figure(3)
    plotcol(q,1./lambda,ind,1,5);
    figure(4)
    plotcol(cos1(q)./lambda,sin1(q)./lambda,ind,1,5);
end

figure(3)
hgexport(gcf,['figures/SPAR_Laplace_' COPULA.name '_GPscale.eps'], formatEPS);
figure(4)
hgexport(gcf,['figures/SPAR_Laplace_' COPULA.name '_limitset.eps'], formatEPS);