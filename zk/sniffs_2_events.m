function [sniff] = sniffs2events(mouse, sess, rec,run,varargin)
%tunrs the sniff into a series of event vectors (sample_stamps)
%input is mouse, sess, rec, run, options
%output is structure with vectors:
%sinff_inh_on   : onset of inhale (zero crossing)
%sniff_exh_on   : onset of exhale (zero crossing)
%sniff_exh_off  : offset of exhale (zero crossing of next)
%sniff_exh_wp   : onset of exhale (parabola estimate, it is the division of the segments for the warped sniffs)


inPar=inputParser;
inPar.addRequired('mouse')
inPar.addRequired('sess')
inPar.addRequired('rec',@ischar)
inPar.addRequired('run',@isscalar)
inPar.parse(mouse,sess, rec,run,varargin{:});
%sniff structures must have been generated with sniff_analysis
%every sniff has
% t0 starting point of the waveform extracted (already inverted if that was
% the case when running sniff_analysis
% 
% t_zer: 3x1 : zero crossings [inh exh inh] that define the two segments of
% the sniff (in ms)
% t_zer_fit: 4x1: zero crossings of the parabola fits [ inhOn inhOff exhOn exhOff] 
% Will keep three numbers that will use for warping:
%sinff_inh_on   : onset of inhale (zero crossing)
%sniff_exh_on   : onset of exhale (zero crossing)
%sniff_exh_off  : offset of exhale (zero crossing of next)
%sniff_exh_wp   : onset of exhale (parabola estimate)

%get the info for the rec and for the run
fn = file_names(mouse, sess,rec);
q = load(fn.sess_info);
recInfo = q.info.rec(strcmpi(rec,{q.info.rec.name}));
runInfo = recInfo.run([recInfo.run.num]==run);

sniff=struct('inh_on',[],'exh_on',[],'exh_off',[],'exh_wp', []);
%get the sniffs for this run
%load the sniff
if ~exist(fn.sniffs,'file')
    warning('Sniffs structure file not found');
    return
end
s=load(fn.sniffs);
%get the sniffs for the run
sniffs=s.sniff([s.sniff.t0]>runInfo.start & [s.sniff.t0]<(runInfo.start+runInfo.duration+1));

for i=1:numel(sniffs)
    t0=sniffs(1).t0;
    sniff(i).inh_on  = sniffs(i).t_zer(1)+t0;
    sniff(i).exh_on  = sniffs(i).t_zer(2)+t0;
    sniff(i).exh_off = sniffs(i).t_zer(3)+t0;
    sniff(i).exh_wp  = sniffs(i).t_zer_fit(2)+t0;
end

sniffsMeta.event='sniff';
sniffsMeta.type ='sniff';

write_m(sniffEvents,sniffsMeta,fn.evt_m,runInfo);


% ~~~~~~~~~~~~~~~~~~~~~~
    function events = write_m(allEvents,allMeta,evFileName,runInfo)
        % open the event file for the record (%yell it if it doesnt
        % exist),
        % save the events and meta under an field 'xx' (two digits for the number of run)
        [evDir,evFile,evExt] = fileparts(evFileName);
        runFieldName     = sprintf('run_%02d',runInfo.num);
        fprintf('\n - Saving sniff into events structure of run %s in matlab file %s%s...\n\t ',runFieldName,evFile,evExt);
        % if the file exists, open it, if not, make it
        evFound = dir(evFileName);
        if ~isempty(evFound)
            fprintf('Found an events file for this rec. I will open it and append the sniff events for this run.\n');
            events = load(evFileName);
            if ~isfield(events,runFieldName)
                warning('Did not find field for events of run %d, not doing anything',runInfo.num);
                return
            end
        else
            fprintf('No previous event file for this rec. Doing nothing.\n');
            return
        end
        
        %add the event structure for this run.
        for iev=1:numel(allMeta)
            meta   = eventsMeta(iev);
            eventRead = allEvents.(meta.event);
            %make every event a horizontal structure with fields
            %eventName_event
            evtFields = fields(eventRead);
            for iF = 1:numel(evtFields);
                fieldName = [meta.event '_' evtFields{iF}];
                events.(runFieldName).(fieldName)           = [eventRead.(evtFields{iF})];
                events.(runFieldName).([fieldName '_meta']) = meta;
            end
        end  
        save(evFileName,'-struct','events');
        fprintf(' - done adding sniff to events file in pr_data folder.\n');
        
    end

end
