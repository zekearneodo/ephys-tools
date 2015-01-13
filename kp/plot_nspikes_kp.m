%% kp 2014-03-04
%% Plots response properties as a function of stimulus parameters
% Currently: just plots for light stimulation, and just number of spikes
%
% Requires that resp structure has been created and saved (if not yet,
% uncomment the call to visualize_responses...
%
% Sessions with and without pulseOffset not resolved yet.


function plot_nspikes_kp(mouse,sess,rec,units,sType)

[~,computerName]=system('hostname');
if ~strcmp(strtrim(computerName),'flipper')
    includePath=fullfile(fileparts(pwd),'current','include');
    addpath(includePath);
end

vr=visualize_responses_kp_10(mouse,sess,rec,units,sType);

unitstr = sprintf('%02d',units);
dm = data_management_tools_032();
fn = dm.file_names(mouse, sess, rec);
fn.resp = fullfile(fn.fold_an_sess,  sprintf('%s_%s_%s_%s_units%s_resp.mat', mouse, sess, rec, sType, unitstr));
load(fn.resp)

% stats_to_plot = {''};

% for istat = 1:length(stats_to_plot)
gs=figure
set(gs,'Nextplot','add')

stims=[resp.stim];

% [stims(i).energy] = [stims.laserDur].*[stims.laserPower];

prop_plot = 'odorConc';
% prop_col_grp = 'laserPower';
% odorplot = [stims.(prop_plot)]*-1
% odors = [stims.(prop_plot)]; latency = [resp.latencyISI]; jitter = resp.jitterISI;
% odors = [odors(3) odors(2) odors(1)]
% latency = [latency(3) latency(2) latency(1)]
% jitter = [jitter(3) jitter(2) jitter(1)]
    
if exist('prop_col_grp')
    stim_grp = unique([stims.(prop_col_grp)])

    grp_colors = {'r' 'g' 'c' 'b' 'm' 'y'};
    hold on
    for igc = 1:numel(stim_grp)
        gInd = [stims.(prop_col_grp)]==stim_grp(igc);
%         errorbar([stims(gInd).(prop_plot)],([resp(gInd).nSpikes_subtBase]),[resp(gInd).nSpkErr],'^','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','k','MarkerEdgeColor',grp_colors{igc})
%         errorbar([stims(gInd).(prop_plot)],([resp(gInd).delay]),[resp(gInd).nSpkErr],'^','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','k','MarkerEdgeColor',grp_colors{igc})
        errorbar([stims(gInd).(prop_plot)],([resp(gInd).latencyISI]),[resp(gInd).jitterISI],'^','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','k','MarkerEdgeColor',grp_colors{igc})
    end
    
else
    disp('No grouping');
    errorbar([-3:-1],([resp.latencyISI]),[resp.jitterISI],'^','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','k')
end
% errorbar([stims.pulseOffset],([resp.nSpk]-[resp.spkBase]),[resp.nSpkErr],'-^','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','k')
xlabel(sprintf('%s',prop_plot));
ylabel('latency');
% xlim([0 120])
% ylim([0 9])

% fnFig=sprintf('/Users/kpenikis/Desktop/%s_%s_%s_%s_units%s.pdf',mouse,sess,rec,sType,unitstr)
% print('-depsc', fnFig);
% end
end

