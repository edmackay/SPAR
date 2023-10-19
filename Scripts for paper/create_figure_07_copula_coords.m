clear
close all

addpath(genpath('..\functions'))
plotsettings_SPAR

%%%%%%%%%%%%%%%%%%%%%% Hua-Joe model linear %%%%%%%%%%%%%%%%%%%%%% 
figure_3panel
xlabel('$u$')
ylabel('$v$')
axis([0 1 0 1])

% lines of equal angle
t=linspace(0,2,100);
for w=0.1:0.1:0.9
    plot(t*(1-w),t*w,'k')
end

% lines of equal radius
w=linspace(0,1,100);
for t=0.2:0.2:1.8
    plot(t*(1-w),t*w,'r')
end

hgexport(gcf, 'figures/Copula_coords_HJ_linear.eps', formatEPS);

%%%%%%%%%%%%%%%%%%%%%% Hua-Joe model log-log %%%%%%%%%%%%%%%%%%%%%% 
figure_3panel
xlabel('$u$')
ylabel('$v$')
set(gca,'xScale','log')
set(gca,'yScale','log')
axis([1e-10 1 1e-10 1])

% lines of equal angle
t=logspace(-20,0,1000);
w0=logspace(-9,-1,9);
for w=[w0 0.5 flip(1-w0)]
    plot(t*(1-w),t*w,'k')
end

% lines of equal radius
w0=logspace(-10,-0.4,1000);
w=[w0 0.5 flip(1-w0)];
for t=logspace(-9,-1,9)
    plot(t*(1-w),t*w,'r')
end

hgexport(gcf, 'figures/Copula_coords_HJ_loglog.eps', formatEPS);

%%%%%%%%%%%%%%%%%%%%%% Wadsworth-Tawn model linear %%%%%%%%%%%%%%%%%%%%%%
figure_3panel
xlabel('$u$')
ylabel('$v$')
axis([0 1 0 1])

% lines of equal angle
t=logspace(-10,0,1000);
for w=0.1:0.1:0.9
    plot(t.^(1-w),t.^w,'k')
end

% lines of equal radius
w=linspace(0,1,1000);
for t=[0.02 0.1:0.1:0.9]
    plot(t.^(1-w),t.^w,'r')
end

hgexport(gcf, 'figures/Copula_coords_WT_linear.eps', formatEPS);

%%%%%%%%%%%%%%%%%%%%%% Wadsworth-Tawn model log-log %%%%%%%%%%%%%%%%%%%%%%
figure_3panel
xlabel('$u$')
ylabel('$v$')
set(gca,'xScale','log')
set(gca,'yScale','log')
axis([1e-10 1 1e-10 1])

% lines of equal angle
t=logspace(-20,0,1000);
for w=0.1:0.1:0.9
    plot(t.^(1-w),t.^w,'k')
end

% lines of equal radius
w=linspace(0,1,1000);
for t=logspace(-18,-2,9)
    plot(t.^(1-w),t.^w,'r')
end

hgexport(gcf, 'figures/Copula_coords_WT_loglog.eps', formatEPS);