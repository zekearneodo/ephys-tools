
function [] = demonstrate_fas_sensitivity(trials,s_test)
sigma = 1;  N = 100; trials = 6; s_test =0;
[f,s,sa_N,ex_cells] = create_population(sigma,N);

for ss = 1:length(s_test)
    
[closest_s, s_index]=min(abs(s-s_test(ss)));
fas_test = f(s_index,:);   % expected response from each neuron to the tested stimulus

figure;
hold on
plot(s,f(:,ex_cells)./50,'LineWidth',3)

r_test = 888*ones(1,N);
for tr = 1:trials
    for nn = 1:N
        r_test(nn,tr,ss) = random('poiss',fas_test(nn));
    end
end


if ss == 1
    plot(sa_N(ex_cells)-0.03,r_test(ex_cells,:,ss)./50,...
        '+k','MarkerSize',16,'LineWidth',2)
else
    plot(sa_N(ex_cells)+0.03,r_test(ex_cells,tr,ss)./50,...
        'ob','MarkerSize',16,'LineWidth',2)
end


plot([s_test(1) s_test(1)],[0 (1.2*max(r_test(:,tr,ss)/50))],':k','LineWidth',2)
if numel(s_test) > 2
    plot([s_test(2) s_test(2)],[0 (1.2*max(r_test(:,tr,ss)/50))],':b','LineWidth',2)
end

xlim([min(s) max(s)])
ylim([0 (1.2*max(r_test(:,tr,ss)/50))])
xlabel('Stimulus','FontSize',28)
ylabel('Normalized response','FontSize',28)
% legend(h2,'S=0','S=0.3')
end
ex_cells;
ex_sa = sa_N(ex_cells);
hold off
end

