function [M,observed_act_vector,P_resp_given_s,log_P_ML] = activity_vector_comp(N,s_test,s_check_prob,lik)
% N = 20;
global rmax s f sa_N

[fas_test,r_test] = present_stimulus(N,s_test);
        observed_act_vector = r_test;
        act_template_true = fas_test;
        max_P_resp_given_s = prod(((act_template_true.^observed_act_vector)./factorial(observed_act_vector)) .* exp(-act_template_true));
        max_log_P_ML = observed_act_vector * log(act_template_true)' * sum(act_template_true);
    
for sch = 1:length(s_check_prob)
    
    [fas_test,r_test] = present_stimulus(N,s_check_prob(sch));
    
    P_resp_given_s(sch) = prod(((fas_test.^observed_act_vector)./factorial(observed_act_vector)) .* exp(-fas_test));
    log_P_ML(sch) = observed_act_vector * log(fas_test)' * sum(fas_test);
    
    figure;
    hold on
    plot(s,f./rmax,'k','LineWidth',0.2);
    plot(sa_N,observed_act_vector./rmax,'.k','MarkerSize',25)
    
    if s_check_prob(sch) == s_test
        [AX,template,prob_r] = plotyy(sa_N,act_template_true./rmax,s_check_prob(1:sch),log_P_ML(1:sch),'plot');
        set(AX(1),'YLim',[0 1.3])
        set(AX(2),'YLim',[0 1.2*max_log_P_ML])
        set(get(AX(1),'Ylabel'),'String','Normalized response');
        set(get(AX(2),'Ylabel'),'String','Probability of observed response given stimulus');
        set(template,'Marker','o','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',10,'LineStyle','none')  
        set(prob_r,'Marker','p','MarkerEdgeColor','b','MarkerSize',10,'MarkerFaceColor','g','LineStyle','--') 
        pause(2)
%         figure;
%         hold on
%         plot(s,f./rmax,'k','LineWidth',0.2);
%         s_test_plot = ones(N).* s_test;
%         plot(s_test_plot,fas_test./rmax,'.k','MarkerSize',25)
%         plot(s_test_plot,observed_act_vector./rmax,'.r','MarkerSize',25)
%         hold off
    else 
%         plot(sa_N,fas_test./rmax,'.r','MarkerSize',25)
        [AX,template,prob_r] = plotyy(sa_N,fas_test./rmax,s_check_prob(1:sch),log_P_ML(1:sch),'plot');
        set(AX(1),'YLim',[0 1.3])
        set(AX(2),'YLim',[0 1.2*max_log_P_ML])
        set(get(AX(1),'Ylabel'),'String','Normalized response');
        set(get(AX(2),'Ylabel'),'String','Probability of observed response given stimulus');
        set(template,'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineStyle','none')  
        set(prob_r,'Marker','p','MarkerEdgeColor','b','MarkerSize',10,'MarkerFaceColor','c','LineStyle','--') 
    end
    M(sch) = getframe;
end

if nargin<4
    hold off
    return
else
    s_check_in = find(ismember(s,s_check_prob));
    P_resp_given_s = ones(length(s_check_prob),1);
    for sch = 1:length(s_check_prob)
        P_resp_given_s(sch) = prod(((f(s_check_in(sch),:).^observed_act_vector)./factorial(observed_act_vector)) .* exp(-f(s_check_in(sch),:)));
    end
disp(P_resp_given_s)

plot(s(s_check_in), P_resp_given_s,'pb','MarkerSize',15,'MarkerFaceColor','c')
    
end
end










