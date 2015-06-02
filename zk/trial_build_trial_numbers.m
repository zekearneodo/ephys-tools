function trial=trial_build_trial_numbers(V,rec,irun)
% build trial structure from a voyeur behavior file already read.
% get the events tables from data_structure_tools->
%V is the voyeur-read file
%trial is the structure of trials
%run is the rec structure to which the trial belongs (e.g: rec(1))
%i is the number of the run
%V is read and trials are created in the trial structure which is returned)
trial = struct();
mouse=V.trial(1).mouse;
sess=V.trial(1).session;
ds = data_structure_tools;
%[rec_info, run_info] = ds.get_info(mouse,sess,rec_name,irun);
s_f = rec.sampling_freq;

%quick and dirty correct case in mouse:
if strcmpi('zkawakem72',mouse)
    mouse = 'ZKawakeM72';
end

%determine whether it is the newer or the older version of the
%Voyeur protocol (the newer version uses Chris's pulsetrain)
if isfield(V.trial,'pulseOnsetDelay_1')
    disp('post-feb2014 protocol detected');
    trial=trial_build(trial);
else
    error('I dont know how to handle protocols without trial numbers')
end

%check if there is a trial correction function to be applied in the
%raw_data folder of the mouse

sess_str = sprintf('%03d',V.trial(1).session);
fn=file_names(mouse,sess,rec.name);
trCFname=fullfile(fn.fold_rd_sess,sprintf('trialCorrect_%s_%s.m',mouse,sess_str));
if exist(trCFname)
    fprintf('Found trial correction function %s\n',trCFname);
    p = path;
    path(p,fn.fold_rd_sess);
    trCorrectFcn=str2func(sprintf('trialCorrect_%s_%s', mouse,sess_str));
    trial=trCorrectFcn(trial);
    path(p);
end

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%
    function tr = trial_build(tr)
        %written for the Feb 2014 protocols of acute1 rig.
        fprintf('Making the trial structures..\n')
        n_tr   = numel(V.trial);
        runStart  = rec.run(irun).start;
        v_events =[];
        if isfield(V,'event')
            events = fields(V.event);
            v_events  = events;
            t_events  = cellfun(@(x) [x 'Times'],events,'UniformOutput',false);
        end
        %get all the events and the good trial numbers for this run
        [trial_numbers, events_list] = ds.make_tables(mouse,sess,rec.name,irun);
        tn_table = trial_numbers.table;
        odor_events = events_list( find(strcmpi('finalValve_1',{events_list.name}))).table.trials;
        laser_events = events_list( find(strcmpi('laser_1',{events_list.name}))).table.trials;
        for kt = 1:numel(tn_table.trialNumber)
            trial_number = tn_table.trialNumber(kt);
            it = find([V.trial.trialNumber]==trial_number);
            if ~(numel(it)==1)
                warning('Wrong number of trials with number %d\n',trial_number);
                continue
            end
            t  = V.trial(it);
            t0_voyeur = t.starttrial;
            t0 = round(tn_table.trialStart(kt)*1000/s_f);
            tr(kt).rec           = rec.name;
            tr(kt).run           = irun;
            tr(kt).runTrialNum   = trial_number;
            tr(kt).runTrialStart = t0;
            tr(kt).recTrialStart = t0+runStart; %trial start from the beginning of the rec (no!!; not all runs begin counting from 0)
            tr(kt).runTrialDur   = round(tn_table.trialEnd(kt)*1000/s_f) - t0;
            
            tr(kt).start      = t0+runStart;  % it's the same as recTrialStart
            tr(kt).duration   = tr(kt).runTrialDur; 
            
            %The timestamps and identity of the stimuli come from the
            %events, but it's a good idea to run double-checks
            
            %%% Odor stimuli
            ot = find([odor_events.trialNumber]==trial_number);
            odor_check = 1;
            conc_check = 1;
            
            %check matching odor and concentration
            if ~isempty(ot)
                odor_check = strcmpi(t.odor,odor_events(ot).odor);
                conc_check = (round(odor_events(ot).odorConcentration*1000)==(round(t.odorconc*1000)));
                
                if odor_check && conc_check
                    tr(kt).odorName   = t.odor;
                    tr(kt).odorConc   = t.odorconc;
                    tr(kt).odorTimes  = round([odor_events(ot).on ; odor_events(ot).off]*1000/s_f) - t0*[1;1];
                else
                    warning('Mismatch in odor stimulus at trialNumber %d',trial_number);
                    tr(kt).odorName   = ['mismatch_' t.odor '_' odor_events(ot).odor];
                    tr(kt).odorConc   = odor_events(ot).odorConcentration;
                    tr(kt).odorTimes=zeros(size(tr(kt).odorTimes));
                end
            else
                tr(kt).odorName   = 'not_an_event';
                tr(kt).odorConc   = t.odorconc;
                tr(kt).odorTimes=[0;0];
            end
            
            if strcmp(t.odor,'dummy') || strcmp(t.odor,'empty') || strcmp(t.odor,'none')
                tr(kt).odorTimes=[0;0];
            end
            %%%
            
            %%% laser stimuli
            lt = find([laser_events.trialNumber]==trial_number);
            %just one pulse!
            % check matching amplitude and Iti
            amplitude_check = 1;
            iti_check       = 1;
            
            if ~isempty(lt)
                iti_check = (t.iti==laser_events(lt).trialIti);
                amplitude_check = (laser_events(lt).laserAmplitude_1==t.amplitude_1);
                
                if iti_check && amplitude_check
                    tr(kt).laserAmp   = t.amplitude_1;
                    tr(kt).laserPower = t.LaserIntensity_1;
                    tr(kt).laserTimes = round([laser_events(lt).on;laser_events(lt).off]*1000/s_f) - t0*[1;1];
                    
                else
                    warning('Mismatch in laser stimulus at trialNumber %d',trial_number);
                    tr(kt).laserAmp   = -1;
                    tr(kt).laserPower = ['mismatch_' t.laserIntensity_1 '_' odor_events(ot).laserIntensity_1];
                    tr(kt).odorTimes=[0;0];
                end
                
            else
                tr(kt).laserAmp   = 0;
                tr(kt).laserPower = [];
                tr(kt).laserTimes = [0;0];
            end
            
            tr(kt).stimID     = t.stimid;
            
            tr(kt).VoyeurParameters = t;
            
            %get the events from the voyeur h5 file, if any
            if t.trialdur>0
                
                t1 = t0_voyeur - 2000;
                t2 = t0_voyeur + 6000;
                for ie = 1:numel(v_events)
                    try
                    in = (V.event.(v_events{ie})(2,:) > t1)&(V.event.(v_events{ie})(1,:) < t2);
                    tr(kt).(t_events{ie}) = V.event.(v_events{ie})(:,in) - t0_voyeur;
                    catch
                    warning('Event %s is screwed up for trial %d, left empty',t_events{ie},kt);
                    tr(kt).(t_events{ie})=[];
                    end
                    
                end
            end
        end
    end %trial_build

end