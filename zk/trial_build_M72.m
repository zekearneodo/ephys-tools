function trial=trial_build_M72(V,rec,irun)
%build trial structure from a voyeur behavior file already read.
%V is the voyeur-read file
%trial is the structure of trials
%run is the rec structure to which the trial belongs (e.g: rec(1))
%i is the number of the run
%V is read and trials are created in the trial structure which is returned)
trial=struct();
%determine whether it is the newer or the older version of the
%Voyeur protocol (the newer version uses Chris's pulsetrain)
if isfield(V.trial,'pulseOnsetDelay_1')
    disp('post-feb2014 protocol detected');
    trial=trial_build_05(trial);
else
    disp('pre-feb2014 protocol detected');
    if isfield(V.trial(1),'mouse') && strcmpi(V.trial(1).mouse,'zkanesthm71') && info.sess==4
        disp('Special case: fvOnTime=0, non-sniff triggered fvalve')
        trial = trial_build_fvonZero(trial);
    else
        trial = trial_build_04(trial);
    end
end

%check if there is a trial correction function to be applied in the
%raw_data folder of the mouse
mouse=V.trial(1).mouse;
sess = sprintf('%03d',V.trial(1).session);
fn=file_names(mouse,sess);
trCFname=fullfile(fn.fold_rd_sess,sprintf('trialCorrect_%s_%s.m',mouse,sess));
if exist(trCFname)
    fprintf('Found trial correction function %s\n',trCFname);
    p = path;
    path(p,fn.fold_rd_sess);
    trCorrectFcn=str2func(sprintf('trialCorrect_%s_%s', mouse,sess));
    trial=trCorrectFcn(trial);
    path(p);
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function tr = trial_build(tr)
        
        n_tr   = numel(V.trial);
        events = fields(V.event);
        start  = rec.run(irun).start;
        
        for it = 1:n_tr
            if (it>1)&&(V.trial(it).trialstart < V.trial(it-1).trialstart)
                continue
            end
            
            kt = it;
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
        runStart  = rec.run(irun).start;
        v_events  = {'sniff_ttl', 'breaks'};
        t_events  = {'sniffTTLTimes', 'breakTimes'};
        laser_on = [1 1 1 1 0 1 1 0 1];
        
        for it = 1:n_tr
            
            kt = it;
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
        runStart  = rec.run(irun).start;
        v_events  = {'sniff_ttl', 'breaks'};
        t_events  = {'sniffTTLTimes', 'breakTimes'};
        laser_on = [1 0 1 1 0];
        %laser_on = [1 1 1 1 1];
        for it = 1:n_tr
            
            kt = it;
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
        runStart  = rec.run(irun).start;
        v_events  = {'sniff_ttl', 'breaks'};
        t_events  = {'sniffTTLTimes', 'breakTimes'};
        laser_on  = zeros(1,numel(rec.run));
        
        for it = 1:n_tr
            
            kt = it;
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
        runStart  = rec.run(irun).start;
        v_events  = {'sniff_ttl', 'breaks'};
        t_events  = {'sniffTTLTimes', 'breakTimes'};
        laser_on  = zeros(1,numel(rec.run));
        for it = 1:n_tr
            kt = it;
            t  = V.trial(it);
            t0 = t.trialstart;
            tr(kt).rec           = rec.name;
            tr(kt).run           = irun;
            tr(kt).runTrialNum   = it;
            tr(kt).runTrialStart = t0;
            tr(kt).recTrialStart = t0+runStart; %trial start from the beginning of the rec (no!!; not all runs begin counting from 0)
            tr(kt).runTrialDur   = t.trialend - t0;
            
            tr(kt).start      = [];             % will be filled after alignment
            tr(kt).duration   = [];
            
            tr(kt).odorName   = t.odor;
            tr(kt).odorConc   = t.odorconc;
            tr(kt).odorTimes  = t.fvOnTime*[1;1] + [0; t.fvdur] - t0*[1;1];
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
        runStart  = rec.run(irun).start;
        v_events  = {'sniff_ttl', 'breaks'};
        t_events  = {'sniffTTLTimes', 'breakTimes'};
        laser_on  = zeros(1,numel(rec.run));
        for it = 1:n_tr
            kt = it;
            t  = V.trial(it);
            t0 = t.trialstart;
            tr(kt).rec           = rec.name;
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
        runStart  = rec.run(irun).start;
        v_events =[];
        if isfield(V,'event')
            events = fields(V.event);
            v_events  = events;
            t_events  = cellfun(@(x) [x 'Times'],events,'UniformOutput',false);
        end
        laser_on  = zeros(1,numel(rec.run));
        for it = 1:n_tr
            kt = it;
            t  = V.trial(it);
            t0 = t.starttrial;
            tr(kt).rec           = rec.name;
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
            end
            
            %just one pulse! (and watch out, voyeur is not saving the
            %number of pulses. Correct that!)
            tr(kt).laserTimes = t.laserontime*[1;1] + [0; t.pulseOnDur_1/1000] - t0*[1;1];
            if t.laserontime==0
                tr(kt).laserTimes=zeros(size(tr(kt).laserTimes));
            end
            tr(kt).laserAmp   = t.amplitude_1;
            tr(kt).laserPower = t.LaserIntensity_1;
            
            tr(kt).stimID     = t.stimid;
            
            tr(kt).VoyeurParameters = t;
            
            if t.trialdur>0
                t1 = t0 - 2000;
                t2 = t0 + 6000;
                for ie = 1:numel(v_events)
                    try
                    in = (V.event.(v_events{ie})(2,:) > t1)&(V.event.(v_events{ie})(1,:) < t2);
                    tr(kt).(t_events{ie}) = V.event.(v_events{ie})(:,in) - t0;
                    catch
                    warning('Event %s is screwed up for trial %d, left empty',t_events{ie},kt);
                    tr(kt).(t_events{ie})=[];
                    end
                    
                end
            end
        end
    end %trial_build_05

end