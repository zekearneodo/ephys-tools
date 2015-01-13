% data management tools:
%   Dima Rinberg
%   Ezequiel Arneodo
%   last updaate:   Feb 13, 2014.
%
%
% collection of the programs for data management in the lab
% version 031 is the same as 03, but instead of using the 'Laser' signal,
% uses the 'TrPin'.
%
% Installation:since it uses the /include of the current
% ephysDataManagement, it's easiest to run it from a folder just one level
% up the local root of the ephysDataManagement.
%trial4 in trial_prep should work fine for awake acute setup experiments
% of jan 2014.
%
% ---------------------------------------------------------------------------------------
% file_names - creates the structure fn for all file names 
% fn = file_names()
% fn = file_names(mouse, sess)
% fn = file_names(mouse_sess,rec);
%
% ---------------------------------------------------------------------------------------
% read_raw_data_info(mouse, sess)
% the program read the following inofmration:
% 1 - log file
% 2 - raw data file directory
% 3 - electrode channel config
% 4 - meta data for individual recordings
% It prints the structure rec for debugging.
% It saves the data into ss_fold, sess_info.mat file
%
% ---------------------------------------------------------------------------------------
% sniffs_to_trial(mouse,sess)
%
% order of post_spike functions:
%  -sniff_analysis_04_zk
%  -afta_da_sorting(outside, have to include here)
%  -trial_prep
%  -trial_shift_estim
%  -trial_alignment 
%  -sniffs_to_trials
%  -spikes_to_trials


function dm = data_management_tools_032()
global dm

[~,computerName]=system('hostname');
if ~strcmp(strtrim(computerName),'flipper')
    includePath=fullfile(fileparts(pwd),'current','include');
    addpath(includePath);
end

    dm.file_names         = @file_names;
    dm.read_raw_data_info = @read_raw_data_info;
    dm.trial_prep         = @trial_prep;
    dm.basicVoyerRead     = @basicVoyerRead;
    dm.trial_shift_estim  = @trial_shift_estim;
    dm.trial_alignment    = @trial_alignment;
    dm.sniffs_to_trials   = @sniffs_to_trials;
    dm.spikes_to_trials   = @spikes_to_trials;
    
    dm.subplot_fig        = @subplot_fig;
    
% pp.read_data_dir  = @read_data_dir;
% pp.ss_prep        = @ss_prep;
% pp.resampling     = @resampling;
end


function trial_prep(mouse, sess)
% reads the voyeur parameters
%builds the trial structure with that data
global trial
    fprintf('trial preparation\n')
    fprintf('============================================================\n')
    fn  = file_names(mouse, sess)
    q   = load(fn.sess_info);
    info = q.info;
    rec  = info.rec;
    %rec=rec(3);

    
    for irec = 1:numel(rec)
        fn = file_names(mouse, sess, rec(irec).name);
        
        trial = struct();
        kt    = 0;
        
        for irun = 1:numel(rec(irec).run)
            fnam = fullfile(fn.fold_rd_sess, rec(irec).run(irun).behav_data);
            fprintf('Rec %s, run %d\n',rec(irec).name,irun);
            V = basicVoyerRead(fnam, 'sniff_ttl','sniff');
            %determine whether it is the newer or the older version of the
            %Voyeur protocol (the newer version uses Chris's pulsetrain)
            if isfield(V.trial,'pulseOnsetDelay_1')
                disp('post-feb2014 protocol detected');
                trial=trial_build_05(trial);
            else
                disp('pre-feb2014 protocol detected');
                if strcmpi(info.mouse,'zkanesthm71') && info.sess==4
                    disp('Special case: fvOnTime=0, non-sniff triggered fvalve')
                    trial = trial_build_fvonZero(trial);
                else
                    trial = trial_build_04(trial);
                end
                
            end
            disp('   Run  trials');
            disp([irun, kt])
        end
        
        save(fn.sess_info, 'info')
        save(fn.trial,     'trial')
    end
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function tr = trial_build(tr)

        n_tr   = numel(V.trial);
        events = fields(V.event);
        start  = rec(irec).run(irun).start;

        for it = 1:n_tr
            if (it>1)&&(V.trial(it).trialstart < V.trial(it-1).trialstart)
                continue
            end

            kt = kt + 1;
            t  = V.trial(it);
            t0 = t.trialstart;
            tr(kt).run       = irun;
            tr(kt).tr_num    = it;
            tr(kt).start     = t0 + start;
            tr(kt).duration  = t.trialend - t0;
            tr(kt).odorTime  = t.fvOnTime*[1;1] + [0; t.fvdur] - t0*[1;1];
            tr(kt).odorName  = t.odor;
            tr(kt).odorConc  = t.odorconc;
            tr(kt).laserTime = t.laserontime*[1;1] + [0; t.duration_1] - t0*[1;1];
            tr(kt).laserAmpl = t.amplitude_1;
            tr(kt).stimID    = t.stimid;

            tr(kt).VoyeurParameters = t;

            if it==n_tr
                t1 = V.event.breaks(2,end);
            else
                t1 = V.trial(it+1).fvOnTime-1;
            end

%             tr(kt).sniffWaveform = V.stream.sniff(t0:t1);

            for ie = 1:numel(events)
                in = (V.event.(events{ie})(2,:) > t0)&(V.event.(events{ie})(1,:)<t1);
                tr(kt).(events{ie}) = V.event.(events{ie})(in) - t0;
            end

        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function tr = trial_build_01(tr)

        n_tr   = numel(V.trial);
%         events = fields(V.event);
        runStart  = rec(irec).run(irun).start;
        v_events  = {'sniff_ttl', 'breaks'};
        t_events  = {'sniffTTLTimes', 'breakTimes'};
        laser_on = [1 1 1 1 0 1 1 0 1];

        for it = 1:n_tr

            kt = kt + 1;
            t  = V.trial(it);
            t0 = t.fvOnTime;
            tr(kt).run           = irun;
            tr(kt).runTrialNum   = it;
            tr(kt).runTrialStart = t0;
            tr(kt).recTrialStart = t0+runStart;
            tr(kt).runTrialDur   = t.trialend - t0;
            
            tr(kt).start      = [];             % will be filled after alignment 
            tr(kt).duration   = [];
            
            tr(kt).odorTimes  = t.fvOnTime*[1;1] + [0; t.fvdur] - t0*[1;1];
            tr(kt).odorName   = t.odor;
            tr(kt).odorConc   = t.odorconc;
            
            tr(kt).laserTimes = t.laserontime*[1;1] + [0; laser_on(irun)*t.duration_1/1000] - t0*[1;1];
            tr(kt).laserAmp   = t.amplitude_1*laser_on(irun);
            tr(kt).laserPower = [];
            
            tr(kt).stimID     = t.stimid;

            tr(kt).VoyeurParameters = t;
            

            if t.trialdur>0
                t1 = t0 - 2000;
                t2 = t0 + 6000;

                for ie = 1:numel(v_events)
                    in = (V.event.(v_events{ie})(2,:) > t1)&(V.event.(v_events{ie})(1,:) < t2);
                    tr(kt).(t_events{ie}) = V.event.(v_events{ie})(:,in) - t0;
                end
            end

        end
    end

    function tr = trial_build_02(tr)
        %written for the october 2013 versions of behavior programs of
        %chronic1 rig.

        n_tr   = numel(V.trial);
        events = fields(V.event);
        runStart  = rec(irec).run(irun).start;
        v_events  = {'sniff_ttl', 'breaks'};
        t_events  = {'sniffTTLTimes', 'breakTimes'};
        laser_on = [1 0 1 1 0];
        %laser_on = [1 1 1 1 1];
        for it = 1:n_tr

            kt = kt + 1;
            t  = V.trial(it);
            t0 = t.trialstart;
            tr(kt).run           = irun;
            tr(kt).runTrialNum   = it;
            tr(kt).runTrialStart = t0;
            tr(kt).recTrialStart = t0+runStart;
            tr(kt).runTrialDur   = t.trialend - t0;
            
            tr(kt).start      = [];             % will be filled after alignment 
            tr(kt).duration   = [];
            
            tr(kt).odorTimes  = t.fvOnTime*[1;1] + [0; t.fvdur] - t0*[1;1];
            tr(kt).odorName   = t.odor;
            tr(kt).odorConc   = t.odorconc;
            
            tr(kt).laserTimes = t.laserontime*[1;1] + [0; laser_on(irun)*t.duration_1/1000] - t0*[1;1];
            tr(kt).laserAmp   = t.amplitude_1*laser_on(irun);
            tr(kt).laserPower = [];
            
            tr(kt).stimID     = t.stimid;

            tr(kt).VoyeurParameters = t;
            

            if t.trialdur>0
                t1 = t0 - 2000;
                t2 = t0 + 6000;

                for ie = 1:numel(v_events)
                    in = (V.event.(v_events{ie})(2,:) > t1)&(V.event.(v_events{ie})(1,:) < t2);
                    tr(kt).(t_events{ie}) = V.event.(v_events{ie})(:,in) - t0;
                end
            end
            
            %Reaconditionning of the trials to:
            %When no odor was presented, odorTimes on and off are equal, and
            %odorConc=0.
            if strcmp(t.odor,'dummy')
                tr(kt).odorTimes(2,:) = tr(kt).odorTimes(1,:);
                tr(kt).odorName   = t.odor;
                tr(kt).odorConc   = t.odorconc;
            end
            
            

        end
    end

    function tr = trial_build_03(tr)
        %written for the december 2013 versions of behavior programs of
        %chronic1 rig.
        %assume no light+odor trials.
        %so if  odor='dummy', laser is on and no odor was present.
        %if odor!='dummy', laser is off and odor was presented
        n_tr   = numel(V.trial);
        events = fields(V.event);
        runStart  = rec(irec).run(irun).start;
        v_events  = {'sniff_ttl', 'breaks'};
        t_events  = {'sniffTTLTimes', 'breakTimes'};
        laser_on  = zeros(1,numel(rec(irec).run));
        
        for it = 1:n_tr
            
            kt = kt + 1;
            t  = V.trial(it);
            t0 = t.trialstart; % start trial from beginning of beh file (run)
            tr(kt).run           = irun; 
            tr(kt).runTrialNum   = it;
            tr(kt).runTrialStart = t0;
            tr(kt).recTrialStart = t0+runStart;
            tr(kt).runTrialDur   = t.trialend - t0;
            
            tr(kt).start      = [];             % will be filled after alignment 
            tr(kt).duration   = [];
            
            tr(kt).odorTimes  = t.fvOnTime*[1;1] + [0; t.fvdur] - t0*[1;1]; 
            tr(kt).odorName   = t.odor;
            tr(kt).odorConc   = t.odorconc;
            
            if strcmp(t.odor,'dummy') || strcmp(t.odor,'empty')
                tr(kt).odorTimes(2,:) = tr(kt).odorTimes(1,:);
                laser_on(irun)=1;
            end

            tr(kt).laserTimes = t.laserontime*[1;1] + [0; laser_on(irun)*t.duration_1/1000] - t0*[1;1];
            tr(kt).laserAmp   = t.amplitude_1*laser_on(irun);
            tr(kt).laserPower = [];
            
            tr(kt).stimID     = t.stimid;
            tr(kt).VoyeurParameters = t;

            if t.trialdur>0
                t1 = t0 - 2000;
                t2 = t0 + 6000;

                for ie = 1:numel(v_events)
                    in = (V.event.(v_events{ie})(2,:) > t1)&(V.event.(v_events{ie})(1,:) < t2);
                    tr(kt).(t_events{ie}) = V.event.(v_events{ie})(:,in) - t0;
                end
            end
        end
    end

    function tr = trial_build_04(tr)
        %written for the jan 2014 protocols of acute1 rig.
        n_tr   = numel(V.trial);
        events = fields(V.event);
        runStart  = rec(irec).run(irun).start;
        v_events  = {'sniff_ttl', 'breaks'};
        t_events  = {'sniffTTLTimes', 'breakTimes'};
        laser_on  = zeros(1,numel(rec(irec).run));
        for it = 1:n_tr
            kt = kt + 1;
            t  = V.trial(it);
            t0 = t.trialstart;
            tr(kt).rec           = rec(irec).name;
            tr(kt).run           = irun;
            tr(kt).runTrialNum   = it;
            tr(kt).runTrialStart = t0;
            tr(kt).recTrialStart = t0+runStart; %trial start from the beginning of the rec (no!!; not all runs begin counting from 0)
            tr(kt).runTrialDur   = t.trialend - t0;
            
            tr(kt).start      = [];             % will be filled after alignment 
            tr(kt).duration   = [];
            
            tr(kt).odorName   = t.odor;
            tr(kt).odorConc   = t.odorconc;
            tr(kt).odorTimes  = t.fvOnTime*[1;1] + [0; 1000] - t0*[1;1]; 
            if strcmp(t.odor,'dummy') || strcmp(t.odor,'empty') || strcmp(t.odor,'none')
                tr(kt).odorTimes(2,:) = tr(kt).odorTimes(1,:);
                laser_on(irun)=1;
            end
                        
            tr(kt).laserTimes = t.laserontime*[1;1] + [0; laser_on(irun)*t.duration_1/1000] - t0*[1;1];
            tr(kt).laserAmp   = t.amplitude_1*laser_on(irun);
            tr(kt).laserPower = [];
            
            tr(kt).stimID     = t.stimid;

            tr(kt).VoyeurParameters = t;
            
            if t.trialdur>0
                t1 = t0 - 2000;
                t2 = t0 + 6000;
                   
                for ie = 1:numel(v_events)
                    if length(V.event.(v_events{ie}))>1
                        in = (V.event.(v_events{ie})(2,:) > t1)&(V.event.(v_events{ie})(1,:) < t2);
                        tr(kt).(t_events{ie}) = V.event.(v_events{ie})(:,in) - t0;
                    end
                end
            end
            
        end
    end %trial_build_04
    

    function tr = trial_build_fvonZero(tr)
        %written for the jan 2014 protocols of acute1 rig.
        n_tr   = numel(V.trial);
        events = fields(V.event);
        runStart  = rec(irec).run(irun).start;
        v_events  = {'sniff_ttl', 'breaks'};
        t_events  = {'sniffTTLTimes', 'breakTimes'};
        laser_on  = zeros(1,numel(rec(irec).run));
        for it = 1:n_tr
            kt = kt + 1;
            t  = V.trial(it);
            t0 = t.trialstart;
            tr(kt).rec           = rec(irec).name;
            tr(kt).run           = irun;
            tr(kt).runTrialNum   = it;
            tr(kt).runTrialStart = t0;
            tr(kt).recTrialStart = t0+runStart; %trial start from the beginning of the rec (no!!; not all runs begin counting from 0)
            tr(kt).runTrialDur   = t.trialend - t0;
            
            tr(kt).start      = [];             % will be filled after alignment 
            tr(kt).duration   = [];
            
            tr(kt).odorName   = t.odor;
            tr(kt).odorConc   = t.odorconc;
            tr(kt).odorTimes  = t.fvOnTime*[1;1] + [0; t.fvdur]; 
            if strcmp(t.odor,'dummy') || strcmp(t.odor,'empty') || strcmp(t.odor,'none')
                tr(kt).odorTimes(2,:) = tr(kt).odorTimes(1,:);
                laser_on(irun)=1;
            end
                        
            tr(kt).laserTimes = t.laserontime*[1;1] + [0; laser_on(irun)*t.duration_1/1000] - t0*[1;1];
            tr(kt).laserAmp   = t.amplitude_1*laser_on(irun);
            tr(kt).laserPower = [];
            
            tr(kt).stimID     = t.stimid;

            tr(kt).VoyeurParameters = t;
            
            if t.trialdur>0
                t1 = t0 - 2000;
                t2 = t0 + 6000;
                   
                for ie = 1:numel(v_events)
                    if length(V.event.(v_events{ie}))>1
                        in = (V.event.(v_events{ie})(2,:) > t1)&(V.event.(v_events{ie})(1,:) < t2);
                        tr(kt).(t_events{ie}) = V.event.(v_events{ie})(:,in) - t0;
                    end
                end
            end
            
        end
    end %trial_build_fvonZero

    function tr = trial_build_05(tr)
        %written for the Feb 2014 protocols of acute1 rig.
        n_tr   = numel(V.trial);
        events = fields(V.event);
        runStart  = rec(irec).run(irun).start;
        v_events  = {'sniff_ttl', 'breaks'};
        t_events  = {'sniffTTLTimes', 'breakTimes'};
        laser_on  = zeros(1,numel(rec(irec).run));
        for it = 1:n_tr
            kt = kt + 1;
            t  = V.trial(it);
            t0 = t.starttrial;
            tr(kt).rec           = rec(irec).name;
            tr(kt).run           = irun;
            tr(kt).runTrialNum   = it;
            tr(kt).runTrialStart = t0;
            tr(kt).recTrialStart = t0+runStart; %trial start from the beginning of the rec (no!!; not all runs begin counting from 0)
            tr(kt).runTrialDur   = t.endtrial - t0;
            
            tr(kt).start      = [];             % will be filled after alignment
            tr(kt).duration   = [];
            
            tr(kt).odorName   = t.odor;
            tr(kt).odorConc   = t.odorconc;
            tr(kt).odorTimes  = t.fvOnTime*[1;1] + [0; t.fvdur] - t0*[1;1];
            if strcmp(t.odor,'dummy') || strcmp(t.odor,'empty') || strcmp(t.odor,'none')
                tr(kt).odorTimes=zeros(size(tr(kt).odorTimes));
            elseif t.fvOnTime == 0
                tr(kt).odorTimes=zeros(size(tr(kt).odorTimes));
            end
            
            %just one pulse! (and watch out, voyeur is not saving the
            %number of pulses. Correct that!)
            tr(kt).laserTimes = t.laserontime*[1;1] + [0; t.pulseOnDur_1/1000] - t0*[1;1];
            if t.laserontime==0
                tr(kt).laserTimes=zeros(size(tr(kt).laserTimes));
            end
            tr(kt).laserAmp   = t.amplitude_1;
            tr(kt).laserPower = str2double(t.LaserIntensity_1(~isletter(t.LaserIntensity_1)));
            tr(kt).pulseOnsetDelay = t.pulseOnsetDelay_1;
            
            tr(kt).stimID     = t.stimid;
            
            tr(kt).VoyeurParameters = t;
            
            if t.trialdur>0
                t1 = t0 - 2000;
                t2 = t0 + 6000;
                
                for ie = 1:numel(v_events)
                    in = (V.event.(v_events{ie})(2,:) > t1)&(V.event.(v_events{ie})(1,:) < t2);
                    tr(kt).(t_events{ie}) = V.event.(v_events{ie})(:,in) - t0;
                end
            end
        end
    end %trial_build_05
    

end

function V = basicVoyerRead(fnam, event_names, stream_names)
%reads the h5 files generated by voyeur.
    if ~iscell(event_names)
        event_names = {event_names};
    end
    if ~iscell(stream_names)
        stream_names = {stream_names};
    end
    
    info = h5info(fnam);
    tr      = read_table(fnam);
    V.trial = tr;
    
  	events = read_events(fnam, event_names);
    for ke = 1:numel(event_names)
        V.event.(event_names{ke}) = events{ke};
    end
    
    [streams, breaks] = read_streams(fnam, stream_names);
    for ks = 1:numel(stream_names)
        V.stream.(stream_names{ks}) = streams(ks,:);
    end
    V.event.breaks = breaks;
    


    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function tr = read_table(fnam)
        
        n_tr = numel(info.Groups);
        table_name = info.Datasets.Name;
        table = h5read(fnam,['/', table_name]);
        ff = fieldnames(table);
        tr = struct();
        
        for kf = 1:numel(ff)
            for it = 1:n_tr
                if ischar(table.(ff{kf}))
                    tr(it).(ff{kf}) = deblank(table.(ff{kf})(:,it)');
                else
                    tr(it).(ff{kf}) = double(table.(ff{kf})(it,1));
                end
            end
        end        
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function [streams, breaks] = read_streams(fnam, sname)
        % read continious sniff waveform and check the packets interuption
        % streams{k} - continious waveforms. if there is a missing packet sn is
        % padded by zeros
        % breaks{k} - [2,nb] - sequence of miisng chnnks of data onset, offset
                
        fprintf('reading all streams (file %s) ...',fnam)
        n_str     = numel(sname);
        sample_on = zeros(1,1e7);
        streams   = zeros(n_str,1e7);
                
        for kt = 1:numel(info.Groups)
            group_name   = info.Groups(kt).Name;
            Events       = h5read(fnam, [group_name,'/Events']);
%             disp(kt)
            for ks = 1:n_str
                streamPackets = h5read(fnam, [group_name,'/',sname{ks}]);
%                 disp(ks)
                for ie = 1:numel(Events.packet_sent_time)
                    nt  = int32(numel(streamPackets{ie}));
%                     disp(ie)
                    ind = Events.packet_sent_time(ie) - nt-1 +(1:nt);
                    if ind == 0
                        continue
                    end
                    streams(ks,ind) = streamPackets{ie};
                    sample_on(ind) = 1;
                end
            end
        end
        
        ind_end = ind(end)+10000;
        ind_end=min(ind_end,numel(streams));
        streams = streams(:,1:ind_end);
        
        breaks_on  = find(diff(sample_on)==-1);
        breaks_off = find(diff(sample_on)==1);
        
        if numel(breaks_on>0)
            %checks if the first on value is greater than the first off,
            %indicating that the sample started on at 1. If the sample started off,
            %then the first on will be the first value in the array.
            if breaks_on(1) > breaks_off(1)
                breaks_on = [1,breaks_on];
            end
            
            % checks if sample is 'on' at the end, if it is, then it forces it
            % off at the last sample.
            if breaks_on(end) > breaks_off(end)
                breaks_off = [breaks_off, ind_end];
            end
        end
       
        breaks = [breaks_on(:)'; breaks_off(:)'];
        fprintf(' done (%d streams) \n',numel(streams))

    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function events = read_events(fnam, ename)
        % reads xx.h5 file trial-by-trail, extract continious sequence ofevents.
        % if the first event ==0, remove this event
        % create a matrix events [2, ne], where
        % ne - number of events
        % event(1,:), event(2,:) are onset and offess of the event
        
        fprintf('reading all events (file %s) ...',fnam)
%         disp(ename);
        % read the initial trial's licks, see if there is an even or odd
        % number of licks. Licks are only transmitted when the beamstate is
        % off. If there is an odd number of lick times, then the first
        % value is off, if even then the first value is on. Want the first
        % value to be on, so discard the first value.
        

        n_names = numel(ename);
        events  = cell(1,n_names);
        for kn = 1:n_names
            events{kn} = zeros(1,1e5);
            ne        = 0;
       
            for kt = 1:numel(info.Groups)
                group_name  = info.Groups(kt).Name;
                eventPacket = h5read(fnam, [group_name,'/', ename{kn}]);
                event_tr    = cell2mat(eventPacket);
                ne_tr       = numel(event_tr);
                events{kn}(ne+(1:ne_tr)) = event_tr;
                ne = ne + ne_tr;
                if (kt == 1)&&mod(ne,2)
                    events{kn} = events{kn}(2:end);
                    ne    = ne-1;
                end
            end
            ne    = floor(ne/2);
            events{kn} = reshape(events{kn}(1:2*ne), 2, ne);
        end
        fprintf(' done (%d events) \n',numel(events{kn}))
    end
    
end

function trial_shift_estim(mouse, sess, rec)
%estimates the misalingment (in trials) of the h5 file of timestapms
%respect to the ephys recordings
%It was made for trigger signals in the beh files that count the time since
%the beginning of the trial.
%For use with the Trial Pin and back compatibility,adapted that the time tt
%does not add any other delay if 'runTrialStart' and 'TrPin' are the
%signals.
%global info

fn = file_names(mouse, sess);
q = load(fn.sess_info);
info = q.info;
if ~isnumeric(sess)
    sess=str2num(sess);
end

fprintf('\n-----------------------------------------------------------\n');
fprintf('Estimating trial shift for mouse %s, session %s...\n',num2str(mouse),num2str(sess));

if nargin<3 || isempty(rec)
    fprintf('\t no record specified, going for all of them...\n');
    reclist = {q.info.rec.name};
else
    reclist = {rec};
end

eTrig_chan = 'TrPin';
bTrig_chan = 'runTrialStart';
trial_shift = 0:10; nts = numel(trial_shift);

for irec = 1:numel(reclist)
    nrec = find(strncmpi(cellstr({info.rec.name}),reclist(irec),1))
    nrun = numel(info.rec(nrec).run);
    fn = file_names(mouse, sess, reclist{irec});
    q = load(fn.trial);
    trial = q.trial;
    q = load(fn.rsm_data, eTrig_chan);

    % beh trial start trig
    if ~ strcmp('runTrialStart',bTrig_chan)
        tt = [1;1]*[trial.runTrialStart] + [trial.(bTrig_chan)];
        bSync = tt(1,:);
    else
        bSync=[trial.runTrialStart];
    end
    
    sf(irec)=figure(), clf
    
    for irun=1:nrun;
        runFig(irun)=figure();
    end
    
    for irun = 1:nrun
        run=info.rec(nrec).run(irun);
        % e-phys trigger
        run_start = info.rec(nrec).run(irun).start;
        run_dur   = info.rec(nrec).run(irun).duration;
        
        % ephys trigger for the run
        eTRun=q.(eTrig_chan)(run.start+(1:run.duration-1));
        if strcmp(mouse,'ZKawakeM72') && sess==1 && nrec==7 && irun==3
            eTRun(1:1.2e5)=0;
        end
        eThresh=min(eTRun)+range(eTRun)*1/2;
        eSyncRunStart = find(diff(eTRun > eThresh)==1);
        %trial endings
        eSyncRunEnd = find(diff(eTRun > eThresh)==-1);
        %Cleanup for non-recorded starts and ends
        if eSyncRunStart(1)>eSyncRunEnd(1)
            eSyncRunEnd(1)=[];
        end
        if eSyncRunStart(end)>eSyncRunEnd(end)
            eSyncRunStart(end)=[];
        end
        eSyncRun=[eSyncRunStart;eSyncRunEnd];
        
        %             figure
        %             plot(eTRun)
        %             plot(diff(eTRun))
        %             hold
        %             plot(eSyncRun,ones(size(eSyncRun))*15000,'ro')
        
        % beh trigger (discard the trials with null duration and zero trialstart)
        in=find( ([trial.run] == irun) & ([trial.runTrialDur]>0) & ([trial.runTrialStart]>0));
        bSyncRun = [bSync(in);bSync(in)+[trial(in).runTrialDur]];
        
        % estimation synchronization error for each trial shift value
        error = zeros(1,nts);
        alpha = zeros(1,nts);
        err_min = 1e10;
        
        for it = 1:nts
            % defining the shifted trigger times
            %positive shift is when the first n trials of the beh were
            %missed

            if trial_shift(it) >= 0
                bTimes = bSyncRun(:,(1+trial_shift(it)):end);
                inShifted= in(1+trial_shift(it):end);
                eTimes = eSyncRun;
            else
                eTimes = eSyncRun;
                bTimes = bSyncRun(:,(1-trial_shift(it)):end);
            end
            
            ne = length(eTimes);
            nb = length(bTimes);
            nn = min([ne, nb]);
            
            %  and alinging triggers by the first trial
            eT0=eTimes(1,1); % the t0 of the aligned ephys, segment, relative to beginning of the run
            bT0=bTimes(1,1);
            eTimes = eTimes(:,1:nn) - eTimes(1,1);
            bTimes = bTimes(:,1:nn) - bTimes(1,1);
            ind    = 1:nn;
            alpha(it) = (eTimes(1,:)/bTimes(1,:));
            error(it) = sqrt(norm(eTimes-alpha(it)*bTimes));
            figure(sf(irec))
            subplot(nrun,3,3*irun)
            plot(trial_shift(it), error(it), 'o', 'Color', 'r', 'LineWidth', 2)
            hold on
            %
            if error(it)<err_min
                err_min = error(it);
                imin    = it;
                figure(sf(irec))
                subplot(nrun,3,3*irun-2)
                plot(eTimes,alpha(it)*bTimes,'x')
                title(trial_shift(it))
                subplot(nrun,3,3*irun-1)
                errorVec=sum(eTimes-alpha(it)*bTimes,1);
                plot(errorVec, 'x')
                figure(runFig(irun))
                clf
                plot(eTimes,ones(size(eTimes))*double(eThresh),'r*')
                hold on
                plot(bTimes*(alpha(it)),ones(size(eTimes))*double(eThresh)*1.2,'mx')
                plot(eTRun(eT0:end))
                title(sprintf('Shift estim + align, rec %s run %02d',info.rec(nrec).name, irun), 'FontSize', 10, 'FontWeight', 'bold')
                trAlign.eTrChan=eTrig_chan;
                trAlign.bTrField=bTrig_chan;
                trAlign.eT0Run = eT0; %start of the aligned ephys rel. to the run.
                trAlign.eTrigsT0 = eTimes; %ephys trig times rel to eT0
                trAlign.iTrialsRun = inShifted(ind); %beh trials that have a matching ephys trial.
            end
            drawnow
        end
        
        trAlign.trialShift          = trial_shift(imin)
        trAlign.frequencyCorrection = alpha(imin);
        info.rec(nrec).run(irun).trAlign=trAlign;
        disp(info.rec(nrec).run(irun).trAlign)
        figure(sf(irec))
        title(sprintf('Shift estim + align, rec %s run %02d',info.rec(nrec).name, irun), 'FontSize', 10, 'FontWeight', 'bold')
    end
end
save(fn.sess_info, 'info')
fprintf('trial shift info saved in %s\n',fn.sess_info)
    
end

function trial_alignment(mouse, sess, rec)
%aligns the trials and the ephys file.
    fn = file_names(mouse, sess);
    q = load(fn.sess_info);
    info = q.info;
    
fprintf('\n-----------------------------------------------------------\n');
fprintf('Aligning trials for mouse %s, session %s...\n',num2str(mouse),num2str(sess));
if nargin<3 || isempty(rec)
    fprintf('\t no record specified, going for all of them...\n');
    reclist = {q.info.rec.name};
else
    reclist = {rec};
end
    
    

    for irec = 1:numel(reclist),    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %use the fields left in the info structure for the run by function trial_shift_estim
        disp(reclist{irec})
        fn = file_names(mouse, sess, reclist{irec});
        
        nrec = find(strncmpi(cellstr({info.rec.name}),reclist(irec),1))
        thisRec=info.rec(nrec);
        disp(thisRec)
        
        % loading behavioral trials for the rec
        q = load(fn.trial);
        trial = q.trial;
        fprintf('rec %s\n ================\n', thisRec.name)
        for irun=1:numel(thisRec.run)
            %load the structure with the alignmet information
            thisRun = thisRec.run(irun)
            disp(thisRun);
            trAlign = thisRun.trAlign;
            
            fprintf('   run %d:\n ________\n', irun)
            
            %by default, set trail.start=0 for all the trials in the run.
            for itr=find([trial.run]==irun)
                trial(itr).start=0;
            end
            %go over all the trials that have an ephys match and relate
            %their start times to those of the ephys run.
            %the indices of the good trials for this run
            iTR=trAlign.iTrialsRun;
            
            for it=1:numel(iTR)
                %there are as many beh trial (numel(iTR)) as ephys triggers
                % shift the beginning of the trial to the corresponding
                % ephys trigger
                %(it)th trigger start relative to the beginning of th rec
                t0=trAlign.eTrigsT0(1,it)+trAlign.eT0Run+thisRun.start - trial(iTR(it)).laserTimes(1)*strcmp(trAlign.eTrChan,'Laser');
                trial(iTR(it)).start=t0;
                fprintf('%d    %d \n', iTR(it), trial(iTR(it)).start)
            end
            
        end %for irun=1:numel(thisRec.run)
        save(fn.trial, 'trial')
        
    end
end

function sniffs_to_trials(mouse, sess, rec)
% the function reads the strucutre sniff.mat from xxx_sniff.mat file and
% sort the sniff parameters into individula trials, save the results tback
% to trial structure to the file xxx_trials.mat
% the mosue and sess must be specifiied. 
% rec is optional
% if rec is specified, it sort sniffs to trials for a specific recors,
% otehrwise, it reads the lsit of all records for a given session
global trial

    fn = file_names(mouse, sess);
    q = load(fn.sess_info);
    
    if nargin<3 || isempty(rec)
        fprintf('\t no record specified, going for all of them...\n');
        rec = {q.info.rec.name};
    else
        rec = {rec};
    end
    disp('tr    num sniffs   delay to 1st inh')
    
    for ir = 1:numel(rec)
        fn = file_names(mouse, sess, rec{ir});
        q = load(fn.trial);
        trial = q.trial; 
        
        q = load(fn.sniffs);
        sniff = q.sniff;
        n_sn  = numel(sniff);
        
        sn_t_zer     = ones(3,1)*[sniff.t0] + [sniff.t_zer];
        
        %ignore if there is no prabola fittings
        if isfield(sniff,'t_zer_fit')
            sn_t_zer_fit = ones(4,1)*[sniff.t0] + reshape([sniff.t_zer_fit], 4, n_sn);
        else
            fprintf('Rec %s has no sniff parabola fittings, leaving sniffParabolaZeroTimes field out\n',rec{ir});
        end

        
        badct = 0;
        for it = 1:numel(trial)
            if isempty(trial(it).start)||(trial(it).start == 0)
                continue
            end
            
            t0 = trial(it).start;
            t1 = t0 - 5000;
            t2 = t0 + 6000;
            
            in_sn = (sn_t_zer(3,:) > t1)&(sn_t_zer(1,:) < t2);
            
            trial(it).sniffZeroTimes      = sn_t_zer([1,2], in_sn) - t0;
            
            if diff([trial(it).odorTimes]) == 0
                continue
            else
                first_sn_in = find(trial(it).sniffZeroTimes(1,:) > trial(it).odorTimes(1),1);
                inh1_delay = trial(it).sniffZeroTimes(1,first_sn_in) - trial(it).odorTimes(1);
                bad_sniff_trials = inh1_delay > 50+diff(trial(it).odorTimes);
                
                if isempty(inh1_delay) || bad_sniff_trials
                    disp('  bad trial: ')
                    trial(it).odorConc = -1; 
                end
            end

            
            if isfield(sniff,'t_zer_fit')
                trial(it).sniffParabZeroTimes = sn_t_zer_fit(:,in_sn)  - t0;
            end
            fprintf('%d          %d          %d \n', it, sum(in_sn),inh1_delay)
        end

        save(fn.trial, 'trial')
    end
end

function spikes_to_trials(mouse, sess, rec)

global trial

    fn = file_names(mouse, sess);
    q = load(fn.sess_info);
    

if nargin<3 || isempty(rec)
    fprintf('\t no record specified, going for all of them...\n');
    rec = {q.info.rec.name};
else
    rec = {rec};
end

    for ir = 1:numel(rec)
        fn = file_names(mouse, sess, rec{ir});
        q = load(fn.trial);
        trial = q.trial;
        
        q = load(fn.spikes);
        
        for iu = 1:numel(q.unit)
            sp{iu} = q.unit(iu).times;
        end
        
        
        for it = 1:numel(trial)
            if isempty(trial(it).start)||(trial(it).start == 0)
                continue
            end
            
            t0 = trial(it).start;
            t1 = t0 - 5000;
            t2 = t0 + 6000;
            
            fprintf('%d    ', it)
            for iu = 1:numel(q.unit)
                trial(it).spikeTimes{iu} = sp{iu}((sp{iu}>t1)&(sp{iu}<t2)) - t0;
                fprintf('%2i   ', numel(trial(it).spikeTimes{iu}));
            end
            
            fprintf('\n')
        end

        save(fn.trial, 'trial');
    end
end

function [gf, gs] = subplot_fig(ny, nx, fig)
    if exist('fig', 'var')
        gf = figure(fig); clf
    else
        gf = figure;
    end
    
    dx = 0.9/nx;    dy = 0.9/ny;
    wx = dx*0.95;   wy = dy*0.95;
    for ix = 1:nx
        for iy = 1:ny
            gs(iy,ix) = subplot('position', [0.08+(ix-1)*dx, 0.95-iy*dy, wx, wy]);
        end
    end
end
