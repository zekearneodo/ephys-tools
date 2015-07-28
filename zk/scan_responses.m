% script to quickly prepare trials and go through cells of a shank making
% histograms.
% Uses packs of functions:
% - post_processing_tools
% - visualize_responses_new
% - trial_build_trial_numbers
% - data_structure_tools

mouse ='ZKawakeM72';
sess = 30;

pp=post_processing_tools();

pp.trial_prep(mouse,sess,'trFunc','trial_build_trial_numbers')
pp.sniffs_to_trials(mouse,sess);
pp.spikes_to_trials(mouse,sess,'NewSpikes',true);

%%
rec = 'e';
%pp.spikes_to_trials(mouse,sess,rec,'NewSpikes',true);
fn = file_names(mouse,sess,rec);

group = 3;
load(fn.spikes);
un = find([unit.group]==group)

close all
for i=1:numel(un)
    vr=visualize_responses_new(mouse,sess,rec,un(i),'odor');
end