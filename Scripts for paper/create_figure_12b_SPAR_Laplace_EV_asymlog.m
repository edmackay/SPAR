clear
close all

plotsettings_SPAR

a=5;
psi1=0.9;
psi2=0.1;
COPULA.name='MEV_asymlog';
COPULA.params=[a,psi1,psi2];

MARGINS.X.name='laplace';
MARGINS.Y.name='laplace';

% calculate origin for polar coordinates
beta=log((1-psi2)/(2*(a-1)*psi1) * (psi1/psi2)^a);
gamma=log((1-psi1)/(2*(a-1)*psi2) * (psi2/psi1)^a);
x0=((a-1)*gamma+a*beta)/(2*a-1);
y0=((a-1)*beta+a*gamma)/(2*a-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate joint distribution
x=linspace(-20,20,500);
[X,Y]=meshgrid(x);
fxy = density_xy(X,Y,COPULA,MARGINS);

% calculate angular density and threshold
q=linspace(-2,2,401);
ur=0*q+10; % set initial threshold at r=10;
u0=ur;
u=0*q;
fq=0*q;
zeta=0*q;
for i=1:length(q)
    disp(num2str(i))
    fq(i)=density_q1(q(i),x0,y0,COPULA,MARGINS);
    if 0<=q(i) && q(i)<=1
        u0(i) = max(abs(x0./cos1(q(i))),abs(y0./sin1(q(i))));
    elseif q(i)>1
        u0(i) = abs(y0/sin1(q(i)));
    elseif -1<=q(i) && q(i)<=0
        u0(i) = abs(x0/cos1(q(i)));
    end
    u(i)=max(u0(i),ur(i));
    if isinf(u(i))
        zeta(i)=0;
    else
        zeta(i)=survivor_rcond_L1(u(i),q(i),x0,y0,COPULA,MARGINS,fq(i));
    end
end

% calculate SPAR approximation
R=abs(X-x0)+abs(Y-y0);
Q=angle_L1(X-x0,Y-y0);
C=abs(cos1(Q));
S=abs(sin1(Q));

FQ=interp1(q,fq,Q,'pchip');
U=interp1(q,u,Q,'pchip');
ZETA=interp1(q,zeta,Q,'pchip');
R(R<U)=NaN;

Q3 = Q>-2 & Q<=-1;
Q4 = Q>-1 & Q<= 0;
Q1 = Q> 0 & Q<= 1;
Q2 = Q> 1 & Q<= 2;

Q1_finger = Q1 & Q>(a-1)/(2*a-1) & Q<a/(2*a-1);
Q1_other = Q1 & ~Q1_finger;

chi=psi2/psi1;
w = S;
z = (1-w).^(a) + (chi*w).^(a);
A = (1-psi1)*(1-w) + (1-psi2)*w + psi1*z.^(1/a);

LAMBDA=0*C;
LAMBDA(Q1_finger)=max(S(Q1_finger),C(Q1_finger))*a - min(S(Q1_finger),C(Q1_finger))*(a - 1);
LAMBDA(Q1_other)=1;
LAMBDA(Q2)=1;
LAMBDA(Q3)=A(Q3);
LAMBDA(Q4)=1;

SIGMA=1./LAMBDA;

fxy_SPAR=ZETA.*FQ.*exp(-(R-U).*LAMBDA).*LAMBDA./(SIGMA+U);

% plot comparison
figure_4panel
contour(x,x,log10(fxy),-60:2:-2)
contour(x,x,log10(fxy_SPAR),-60:2:-2,'k--')
plot(x0+u.*cos1(q),y0+u.*sin1(q),'r','LineWidth',2)
plot([x0 x0],[-20 20],'k-')
plot([-20 20],[y0 y0],'k-')
axis([-1 1 -1 1]*20)
scatter(x0,y0,'ro','filled')
set(gca,'YTick',-20:10:20)
set(gca,'xTick',-20:10:20)
title(['$\alpha=' num2str(a) '$, $\gamma_1=' num2str(psi1) '$, $\gamma_2=' num2str(psi2) '$'] )
xlabel('$x$')
ylabel('$y$')
hgexport(gcf,['figures/SPAR_Laplace_' COPULA.name '_XY_density.eps'], formatEPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate angular density
figure_4panel
xlabel('$q$')
ylabel('$f_Q(q)$')

q=linspace(-2,2,401);
ind=0;
for a=[1.5 2 3 5]
    ind=ind+1;
    COPULA.params=[a,psi1,psi2]; 
    fq=0*q;
    for i=1:length(q)
        disp(num2str(i))
        fq(i)=density_q1(q(i),x0,y0,COPULA,MARGINS);
    end
    h=plotcol(q,fq,ind,1,5);
    h.DisplayName=['$\alpha=' num2str(a) '$'];
end
L=legend;
L.Box='off';
L.Location="northwest";
hgexport(gcf,['figures/SPAR_Laplace_' COPULA.name '_angular_density.eps'], formatEPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate limit sets and rate parameter
figure_4panel
axis([-2 2 0 2])
xlabel('$q$')
ylabel('$\sigma(q)$')

figure_4panel
axis([-1 1 -1 1])
xticks(-1:1)
yticks(-1:1)
xlabel('$x$')
ylabel('$y$')

q=linspace(-2,2,401);
q3 = q>=-2 & q<=-1;
q4 = q>-1 & q<= 0;
q1 = q> 0 & q<= 1;
q2 = q> 1 & q<= 2;
c=abs(cos1(q));
s=abs(sin1(q));

ind=0;
for a=[1.5 2 3 5]
    ind=ind+1;

    lambda=0*q;
    if a==1
        lambda=1+0*q;
    else
        q1_finger = q1 & q>(a-1)/(2*a-1) & q<a/(2*a-1);
        q1_other = q1 & ~q1_finger;

        chi=psi2/psi1;
        w = s;
        z = (1-w).^a + (chi*w).^a;
        A = (1-psi1)*(1-w) + (1-psi2)*w + psi1*z.^(1/a);

        lambda(q1_finger)=max(s(q1_finger),c(q1_finger))*a - min(s(q1_finger),c(q1_finger))*(a - 1);
        lambda(q1_other)=1;
        lambda(q2)=1;
        lambda(q3)=A(q3);
        lambda(q4)=1;
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