%% troubleshooting things

addpath /Users/kpenikis/Documents/MATLAB/ephysDataManagement/kp
addpath /Users/kpenikis/Documents/Rinberg/Experiment
% 
% dm = data_management_tools_032_kp

mouse = 'KPawakeM72';
sess = '003';

load /Users/kpenikis/Documents/Rinberg/Experiment/pr_data/mKPawakeM72_003/KPawakeM72_003_c_trial.mat
load /Users/kpenikis/Documents/Rinberg/Experiment/pr_data/mKPawakeM72_003/KPawakeM72_003_c_rsm.mat
load /Users/kpenikis/Documents/Rinberg/Experiment/pr_data/mKPawakeM72_003/KPawakeM72_003_c_sniff.mat
Sniff_inv = Sniff';

for ii = 1:180
    odorName{ii} = trial(ii).odorName;
end

et = strfind(odorName,'ethyl_tiglate');

jj = 1;
for ii = 1:180
    if et{ii} == 1
        et_in(jj) = ii;
        jj = jj+1;
    end
end

et_trialstart = [trial(et_in).start];

for ii = 1:length(et_in)
    et_FV(ii) = et_trialstart(ii) + trial(et_in(ii)).odorTimes(1);    
    trig_point(ii) = Sniff_inv(et_FV(ii));
    
end
trig_point = double(trig_point);
sd_trig = std(trig_point);
med_trig = median(trig_point);
bad_trig = find(trig_point<med_trig-sd_trig);


% figure;
% hist(trig_point)

pick_tr = 21;

sniff_zeros = [trial(et_in(pick_tr)).sniffZeroTimes(1,:) trial(et_in(pick_tr)).sniffZeroTimes(2,:)];
sniff_zeros = sniff_zeros + et_trialstart(pick_tr);
window = 3000;

figure;
hold on
plot((et_FV(pick_tr)-window) : (et_FV(pick_tr)+window), Sniff_inv((et_FV(pick_tr)-window) : (et_FV(pick_tr)+window)))
plot([et_FV(pick_tr) et_FV(pick_tr)], [-1e4 1e4],'--g','MarkerSize',25)
plot(sniff_zeros,zeros(length(sniff_zeros)),'.k')
xlim([et_trialstart(pick_tr)-window et_trialstart(pick_tr)+window])
plot(et_FV(pick_tr),trig_point(pick_tr),'sr')


figure;
hold on
bar(trig_point)
plot([0 length(trig_point)],[-sd_trig -sd_trig],'--g')


