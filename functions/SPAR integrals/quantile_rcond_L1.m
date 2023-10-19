function r=quantile_rcond_L1(Q,q1,x0,y0,COPULA,MARGINS,fq1)

opts=optimset('display','none');
r=fsolve(@fun,1,opts);

% R=linspace(5,15,50);
% for i=1:length(R)
%     dif(i)=fun(R(i));
% end
% figure
% plot(R,dif)
% set(gca,'yscale','log')


    function dif=fun(R)
        dif = abs(Q - survivor_rcond_L1(R,q1,x0,y0,COPULA,MARGINS,fq1))/Q;
    end
end