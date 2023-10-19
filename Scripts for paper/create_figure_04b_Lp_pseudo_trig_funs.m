clear
close all

plotsettings_SPAR

figure('Position',[100 100 766 420])
tiledlayout(2,1)

ax1=nexttile;
hold on; grid on
xlabel('$q$')
ylabel('$\cos_p(q)$')
set(gca,'FontSize',14)

ax2=nexttile;
hold on; grid on
xlabel('$q$')
ylabel('$\sin_p(q)$')
set(gca,'FontSize',14)

p=[0.5 1 2 4];
colmax=length(p)+2;
for i=1:length(p)

    % calculate circumference of Lp unit circle
    w0=2^(-1/p(i));
    Cp=8*integral(@(w)integrand_fun(w,p(i)),0,w0);

    % calculate q for w0<=x<=1
    wx=linspace(1,w0,100);
    wy=(1-wx.^p(i)).^(1/p(i));
    q=0*wx;
    for j=2:length(wx)
        q(j)=integral(@(w)integrand_fun(w,p(i)),0,wy(j))*(4/Cp);
    end
        
    % plot pseudo-trig functions using symmetry
    nexttile(1)
    plotcol(q,wx,i,1,colmax,'-')       % Q1
    plotcol(1-q,wy,i,1,colmax,'-')     % Q1
    plotcol(-q,wx,i,1,colmax,'-')      % Q4
    plotcol(q-1,wy,i,1,colmax,'-')     % Q4
    plotcol(2-q,-wx,i,1,colmax,'-')    % Q2
    plotcol(1+q,-wy,i,1,colmax,'-')    % Q2
    plotcol(q-2,-wx,i,1,colmax,'-')    % Q3
    plotcol(-1-q,-wy,i,1,colmax,'-')   % Q3

    nexttile(2)
    h=plotcol(q,wy,i,1,colmax,'-');     h.DisplayName=['$p=' num2str(p(i)) '$'];
    h=plotcol(1-q,wx,i,1,colmax,'-');   h.HandleVisibility='off';
    h=plotcol(-q,-wy,i,1,colmax,'-');   h.HandleVisibility='off';
    h=plotcol(q-1,-wx,i,1,colmax,'-');  h.HandleVisibility='off';
    h=plotcol(2-q,wy,i,1,colmax,'-');   h.HandleVisibility='off';
    h=plotcol(1+q,wx,i,1,colmax,'-');   h.HandleVisibility='off';
    h=plotcol(q-2,-wy,i,1,colmax,'-');  h.HandleVisibility='off';
    h=plotcol(-1-q,-wx,i,1,colmax,'-'); h.HandleVisibility='off';
end

% manually plot L^infinity trig functions
nexttile(1)
plotcol([-2 -1.5],[-1 -1],i+1,1,colmax,'-')
plotcol([-1.5 -0.5],[-1 1],i+1,1,colmax,'-')
plotcol([-0.5 0.5],[1 1],i+1,1,colmax,'-')
plotcol([0.5 1.5],[1 -1],i+1,1,colmax,'-')
plotcol([1.5 2],[-1 -1],i+1,1,colmax,'-')

nexttile(2)
h=plotcol([-2 -1.5],[0 -1],i+1,1,colmax,'-');    h.DisplayName='$p=\infty$';
h=plotcol([-1.5 -0.5],[-1 -1],i+1,1,colmax,'-'); h.HandleVisibility='off';
h=plotcol([-0.5 0.5],[-1 1],i+1,1,colmax,'-');   h.HandleVisibility='off';
h=plotcol([0.5 1.5],[1 1],i+1,1,colmax,'-');     h.HandleVisibility='off';
h=plotcol([1.5 2],[1 0],i+1,1,colmax,'-');       h.HandleVisibility='off';

L=legend;
L.Box='off';
L.Layout.Tile = 'East';

hgexport(gcf, 'figures/Pseudotrig_funs.eps', formatEPS);


function val=integrand_fun(v,p)
    vp=v.^p;
    term=(vp./(1-vp)).^(2-2/p);
    val=sqrt(1+term);
end