
function [ML,MAP] = estimate_error(s_test,trials)
% s_test = -5:0.25:5;
% trials = 60;
% N = 100;
% sigma = 1;
% sML = zeros(length(s_test),trials);
sMAP = zeros(length(s_test),trials);
parfor sti = 1:length(s_test)
    
    si_test = s_test(sti);
    
%     sML(sti,:) = ML_estimate(sigma,N,si_test,trials);
    
    sMAP(sti,:) = MAP_estimate(sigma,N,si_test,trials);

end

for sti = 1:length(s_test)
% bias of the models
si_test = s_test(sti);
mean_ML = mean(sML(sti,:));
mean_MAP = mean(sMAP(sti,:));

ML.b_est(sti) = mean_ML - si_test;
MAP.b_est(sti) = mean_MAP - si_test;

% variance of estimator
ML.var_est(sti) = mean((sML(sti,:) - mean_ML).^2);
MAP.var_est(sti) = mean((sMAP(sti,:) - mean_MAP).^2);

% squared estimation error
ML.sq_err(sti) = ML.var_est(sti) + ML.b_est(sti);
MAP.sq_err(sti) = MAP.var_est(sti) + MAP.b_est(sti);
end

end