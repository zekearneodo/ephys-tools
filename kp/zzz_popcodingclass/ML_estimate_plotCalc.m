
function [sML,fa,ra,P_ra_s,sa_N_view] = ML_estimate_plotCalc(sigma,N,s_test,trials)
T=1; figure; hold on
view_neuron_N = [24 38 48];
for vn = 1:length(view_neuron_N)
    view_N = view_neuron_N(vn);
for tr = 1:trials
    [r_test,f,s,sa_N] = present_stimulus(sigma,N,s_test);
    sa_N_view(vn) = sa_N(view_N);
    
    for ii = 1:length(s)
        log_P_ML(ii,:) = T * sum(r_test .* log(f(ii,:)));
    end
    
        if tr < 8'
            if tr < 2
                
                if vn==1
            plot(s,log_P_ML./max(log_P_ML),'g','LineWidth',3)
                end
            
            plot(sa_N(view_N),r_test(view_N)./50,'.k','MarkerSize',30) % observed ra(s)
            plot(sa_N(view_N),f(51,view_N)./50,'.m','MarkerSize',30)  % predicted fa(s)
            plot(s,f(:,view_N)./50,':m','LineWidth',2)  % tuning curve of neuron a
            xlim([min(s) max(s)])
            ylim([-.25 1.2*max(r_test)/50])
            plot([s_test(1) s_test(1)],[-.25 1.2*max(r_test)/50],'--k')
            end
    
        fa(tr,vn) = f(51,view_N);
        ra(tr,vn) = r_test(view_N);
        r_ratio(tr,vn) = ra(tr,vn)/fa(tr,vn);
        P_ra_s(tr,vn) = ((fa(tr,vn)^ra(tr,vn))/factorial(ra(tr,vn))) * exp(-fa(tr,vn));
        
        end
        
    sML(tr) = s(log_P_ML==max(log_P_ML));

    %     sML_direct(tr) = sum(r_test.*sa_N) / sum(r_test);
end
    mean_p = mean(P_ra_s(:,vn))
%     scatter(sa_N(view_N),mean(P_ra_s(vn,:)),20,'y','p','fill')
    scatter(sa_N(view_N),mean_p,100,'y','fill')
    

end
xlabel('Stimulus','FontSize',28)
ylabel('P(r|s)','FontSize',28)
end


