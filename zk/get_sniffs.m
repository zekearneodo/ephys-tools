
load('/experiment/export_data/ZKawakeM72_013_e_sniff.mat')
load('/experiment/export_data/ZKawakeM72_013_e_noStimSniff.mat')

i_s = 1;


sn=noStimSniffs(is)
figure
plot(sn.waveform)
hold on
t0 = sn.t0;
plot(-Sniff(sn.t0:sn.t0+numel(sn.waveform)),'r.')