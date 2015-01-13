%general function for correcting trial structure
%change odor names and vial concentrations for mouse KPawakeM72, sess 003.
%gets a trial structure, applies the correction, and retunrs it.
%it is meant to be placed in the pr_data folder.
%when the first function that prepares the trial structure finds it in the
%folder, it runs it.
%the function can be whatever, as long as it doesn't change field names of
%the trial structure.
%Careful, make it so that if you run the function several times it won't
%iterate a correction to screw things up!

function trial=trial_correct_KPawakeM72(trial)
%here all the corrections that apply to all the trials in the sess.

fprintf('Correcting trials for session 003 of KPanesthM72...\n');
fprintf('\t(correct odor names and recalculate delivered concentrations,\n necessary due to switch in voyeur_rig_config files');

vial(4).odor='pinene'
vial(4).vialConc=0.01
vial(5).odor='4-methylacetophenone';
vial(5).vialConc=0.01;
vial(6).odor='methyl_salicylate';
vial(6).vialConc=0.01;
vial(7).odor='1-octanol';
vial(7).vialConc=0.01;

%go through all the vials making the corrections in the trials that involve
%those vials.

for iv=1:numel(vial)
    if isempty(vial(iv).name)
        continue
    end
    
    %get the trials
    vp=[trial.VoyeurParameters];
    ivTrials=find([vp.vial]==iv);
        
    %apply correction to every trial selected
    for it=ivTrials
        trial(it).odorName=vial(iv).odor;
        trial(it).odorConc=vp(it).nitrogen/(vp(it).nitrogen+vp(it).air)/vp(it).dillution*vial(iv).vialConc;
    end
end

end