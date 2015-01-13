
function [sML,sMAP] = MAP_estimate(sigma,N,s_test,trials)
T=1;rmax=50;

[sML] = ML_estimate(sigma,N,s_test,trials)

for tr = 1:trials
    [r_test,f,s,sa_N] = present_stimulus(sigma,N,s_test);
    Ps = pdf('Normal',s,-2,1);
    
    for ii = 1:length(s)
        log_P_MAP(ii,:) = T * sum(r_test .* log(f(ii,:))) + log(Ps(ii));
    end

    sMAP(tr) = s(log_P_MAP==max(log_P_MAP));





[closest, sa_index]=min(abs(sa_N-s_test));
exTC = f(:,sa_index);
if tr ==1
%     plot(s,Ps,'g')
    plot(s,exTC./rmax,'--k')
%     plot(s,log_P_MAP./max(log_P_MAP),'b','LineWidth',2)
%     scatter(sa_N, r_test./rmax,'k')
%     ylim([0 1.2])
tit = sprintf('N neurons = %i, sigma = %i, trials = %i',N,sigma,trials);
title(tit)
end
end
% aaa=3;
end
