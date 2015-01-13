
function [log_P_ML,sML] = ML_estimate(N,s_test,trials)

global s rmax f sa_N
for tr = 1:trials
    [fas_test,r_test] = present_stimulus(N,s_test);
    
    for ii = 1:length(s)
        log_P_ML(ii,:) = r_test * log(f(ii,:))' * sum(f(ii,:));
    end
    
    sML(tr) = s(log_P_ML==max(log_P_ML));
    sML_direct(tr) = sum(r_test.*sa_N) / sum(r_test);
    
    avg_log_P_ML = mean(log_P_ML,2);
    
    figure;
    hold on
    plot(s,log_P_ML./max(log_P_ML),'m','LineWidth',2)
    plot(sa_N(1:N),r_test./rmax,'ok')
    xlim([min(s) max(s)])
    ylim([-0.5 1.5])


end
end


