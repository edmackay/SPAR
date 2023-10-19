clear

plotsettings_SPAR

origin=2; %%%%%%%%% set to 1 or 2 to show effect of different origins

% copula parameters
a=5;
psi1=0.9;
psi2=0.1;

if origin==1
    x0=0;
    y0=0;
elseif origin==2
    % Choice of origin that ensures T2=T3 at fixed angle
    beta=log((1-psi2)/(2*(a-1)*psi1) * (psi1/psi2)^a);
    gamma=log((1-psi1)/(2*(a-1)*psi2) * (psi2/psi1)^a);
    x0=((a-1)*gamma+a*beta)/(2*a-1);
    y0=((a-1)*beta+a*gamma)/(2*a-1);
end

figure_3panel
set(gca,'yscale','log')
set(gca,'yminortick','off')
set(gca,'xtick',0.3:0.2:0.8)
xlabel('$q$')
ylabel('$\delta_1(r,q,c)$')
title(['$(x_0,y_0)=(' num2str(x0,3) ',' num2str(y0,3) ')$'])
xlim([0.2,0.8])

w=linspace(0.2,0.8,1000);

ind=0;
for r=10:5:30
    ind=ind+1;
    
    u = margin_cdf_laplace(x0+r*w);
    v = margin_cdf_laplace(y0+r*(1-w));
    c = copula_density_MEV_asymlog(u,v,[a,psi1,psi2]);
    
    h=plotcol(w,c,ind,1,6);
    h.DisplayName=['$r=' num2str(r) '$'];
    
end
if origin==2
    L=legend;
    L.Box='off';
end
ylim([1e-2 1e7])
set(gca,'YTick',10.^(-2:2:6))
ax=gca;
ax.Children=ax.Children([5 4 3 2 1]);

hgexport(gcf, ['figures/AsymLog_origin_' num2str(origin) '.eps'], formatEPS);