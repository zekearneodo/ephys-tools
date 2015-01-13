function plot_nspikes()

dm = data_management_tools_031();
fn = dm.file_names(mouse, sess, rec, stat);

%load the resp file
load('/home/zeke/data/analysis/mouse_KPawakeM72/sess_004/KPawakeM72_004_a_laser_units03091011121314_resp.mat');
fnFig='/home/zeke/Dropbox/senseofsmell/drafts/r01/delaySpikes/spikesVsDelay.eps';

gs=figure
set(gs,'Nextplot','add')

stims=[resp.stim];

errorbar([stims.pulseOffset],([resp.nSpk]-[resp.spkBase]),[resp.nSpkErr],'-^','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','k')
xlabel('delay (ms)');
ylabel('spkes evoked');
xlim([0 120])
ylim([0 9])

print('-depsc2',fnFig);

end

