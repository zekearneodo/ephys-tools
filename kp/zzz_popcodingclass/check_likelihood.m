
function [P_ra_s,sML] = check_likelihood(N,s_check_prob,fas)

global s rmax f sa_N
    
    for ii = 1:length(s_check_prob)
        P_ra_s(ii) = ((fa(tr,vn)^ra(tr,vn))/factorial(ra(tr,vn))) * exp(-fa(tr,vn));
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


