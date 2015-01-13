
function [fas_test,r_test] = present_stimulus(N,s_test)

global s sa_N f
[sa_N,f] = create_population(N);

for ss = 1:length(s_test)
    
[closest_s, s_index]=min(abs(s-s_test(ss)));
fas_test = f(s_index,:);   % expected response from each neuron to the tested stimulus

r_test = 666*ones(1,N);
for nn = 1:N
    % spiking response from each neuron in population
    % to the presented stimulus
    r_test(nn) = random('poiss',fas_test(nn));
end

% if ss == 1
%     plot(sa_N(ex_cells)-0.075,r_test(ex_cells)./50,'ok')
% else
%     plot(sa_N(ex_cells)+0.075,r_test(ex_cells)./50,'ob')
% end
end

% ylim([0 1.2*max(r_test)/50])
% plot([s_test(1) s_test(1)],[0 1.2*max(r_test)/50],'--k')
end

