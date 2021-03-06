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


function pp = post_processing_tools()
%global pp

[~,computerName]=system('hostname');
if ~strcmp(strtrim(computerName),'flipper')
    includePath=fullfile(fileparts(pwd),'current','include');
    addpath(includePath);
end

    %pp.file_names         = @file_names;
    %file names is in the /include which is in the path of all users of the
    %CS
     
    pp.trial_prep         = @trial_prep;
    pp.basicVoyerRead     = @basicVoyerRead;
    pp.trial_shift_estim  = @trial_shift_estim;
    pp.trial_alignment    = @trial_alignment;
    pp.sniffs_to_trials   = @sniffs_to_trials;
    pp.spikes_to_trials   = @spikes_to_trials;
    pp.trials_match       = @trials_match;
    pp.match_trials       = @match_trials;
    pp.subplot_fig        = @subplot_fig;
    
end


function trial_prep(mouse, sess,varargin)
% trial_prep(mouse,sess,'trFunc',function name,'trFuncPath',function path)
% reads the voyeur parameters
%builds the trial structure with that data
%the inputs:
%trFuncName:    function to build the trial structure after reading the voyeur
%               parameters (default is 'trial_build_M72')
%trFuncPath:    path of the function (default is
%               '/experiment/ephysDataManagement/zk/)

inPar=inputParser;
addRequired(inPar,'mouse')
addRequired(inPar,'sess')



addParameter(inPar,'trFunc','trial_build_M72',@ischar)
addParameter(inPar,'trFuncPath',fullfile(get_local_disk,'ephysDataManagement','zk'),@ischar)
addParameter(inPar, 'recList', {}, @iscell)

parse(inPar,mouse,sess,varargin{:})

%for debugging
global trial

addpath(inPar.Results.trFuncPath);
trialBuildFunction=str2func(inPar.Results.trFunc);

fprintf('trial preparation\n')
fn  = file_names(mouse, sess)
fprintf('MOUSE %s - SESS %s\n',fn.mouse,fn.sess)
fprintf('============================================================\n')
q   = load(fn.sess_info);

info = q.info;

if isempty(inPar.Results.recList)
    rec  = info.rec;
else
    rec = cellfun(@(x) info.rec(strcmpi(x,{info.rec.name})), inPar.Results.recList);
end

%rec=rec(3);


for irec = 1:numel(rec)
    fn = file_names(mouse, sess, rec(irec).name);
    
    trial = []; %starts as an empty memory space, to be filled by subsequent calls of trFunc
        
    for irun = 1:numel(rec(irec).run)
        thisRun=rec(irec).run(irun);
        fnam = fullfile(fn.fold_rd_sess, thisRun.behav_data);
        fprintf('Rec %s, run %d\n',rec(irec).name,irun);
        V = basicVoyeurRead(fnam, {'sniff_ttl' 'lick1' 'lick2'},'sniff');
        trial=[trial trialBuildFunction(V,rec(irec),irun)]; %#ok<AGROW>
        disp('   Run  trials');
        disp([irun, numel(trial)])
    end
    
    save(fn.sess_info, 'info')
    save(fn.trial,     'trial')
end
    
end

function V = basicVoyeurRead(fnam, event_names, stream_names)
%reads the h5 files generated by voyeur.
    if ~iscell(event_names)
        event_names = {event_names};
    end
    if ~iscell(stream_names)
        stream_names = {stream_names};
    end
    
    hinfo = h5info(fnam);
    tr      = read_table(fnam);
    V.trial = tr;
    
    %remove the events that are attempted to be read but are not actually
    %in the file
    presentEvents=cellfun(@(x) any(strcmpi(x,{hinfo.Groups(1).Datasets.Name})),event_names)
    if numel(presentEvents>0)
        warning('There were %d events not found in the h5 file',numel(presentEvents))
        event_names=event_names(find(presentEvents));
    end
    
    if ~isempty(event_names)
        events = read_events(fnam, event_names);
        for ke = 1:numel(event_names)
            V.event.(event_names{ke}) = events{ke};
        end
    end
    [streams, breaks] = read_streams(fnam, stream_names);
    for ks = 1:numel(stream_names)
        V.stream.(stream_names{ks}) = streams(ks,:);
    end
    V.event.breaks = breaks;
    


    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function tr = read_table(fnam)
        
        n_tr = numel(hinfo.Groups);
        table_name = hinfo.Datasets.Name;
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
                
        for kt = 1:numel(hinfo.Groups)
            group_name   = hinfo.Groups(kt).Name;
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
       
            for kt = 1:numel(hinfo.Groups)
                group_name  = hinfo.Groups(kt).Name;
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

function trial_shift_estim(mouse, sess, varargin)
%estimates the misalingment (in trials) of the h5 file of timestapms
%respect to the ephys recordings
%It was made for trigger signals in the beh files that count the time since
%the beginning of the trial.
%For use with the Trial Pin and back compatibility,adapted that the time tt
%does not add any other delay if 'runTrialStart' and 'TrPin' are the
%signals.
%global info


inPar=inputParser;
inPar.addRequired('mouse')
inPar.addRequired('sess')
inPar.addParamValue('figures','plot',@(x) ischar(x) && any(strcmpi(x,{'plot','noplot'})));
inPar.addOptional('rec','',@(x)iscell(x) || (ischar(x) && ~any(strcmpi(x,inPar.Parameters))) );

inPar.parse(mouse,sess, varargin{:});
figures=inPar.Results.figures;

fn = file_names(mouse, sess);
q = load(fn.sess_info);
info = q.info;
if ~isnumeric(sess)
    sess=str2num(sess);
end

fprintf('\n-----------------------------------------------------------\n');
fprintf('Estimating trial shift for mouse %s, session %s...\n',num2str(mouse),num2str(sess));

rec=inPar.Results.rec;
if nargin<3 || isempty(rec)
    fprintf('\t no record specified, going for all of them...\n');
    reclist = {q.info.rec.name};
else
    if iscell(rec)
        reclist =rec;
    else
        reclist={rec}
    end
end

eTrig_chan = 'TrPin';
bTrig_chan = 'runTrialStart';
trial_shift = 0:4; nts = numel(trial_shift);

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
    
    
    
    if strcmpi(figures,'plot')
        sf(irec)=figure(), clf
        for irun=1:nrun;
            runFig(irun)=figure();
        end
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
        
%                     figure
%                     plot(eTRun)
%                     plot(diff(eTRun))
%                     hold
%                     plot(eSyncRun,ones(size(eSyncRun))*15000,'ro')
        
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
            if strcmpi(figures,'plot')
                figure(sf(irec))
                subplot(nrun,3,3*irun)
                plot(trial_shift(it), error(it), 'o', 'Color', 'r', 'LineWidth', 2)
                xlabel('shift');
                ylabel('error sum');
                hold on
            end
            %
            if error(it)<err_min
                err_min = error(it);
                imin    = it;
                if strcmpi(figures,'plot')
                    figure(sf(irec))
                    subplot(nrun,3,3*irun-2)
                    plot(eTimes,alpha(it)*bTimes,'x')
                    xlabel('ephys trig');
                    ylabel('behav trig');
                    title(trial_shift(it))
                    subplot(nrun,3,3*irun-1)
                    errorVec=sum(eTimes-alpha(it)*bTimes,1);
                    plot(errorVec, 'x')
                    xlabel('trial');
                    ylabel('error');
                    figure(runFig(irun))
                    clf
                    plot(eTimes,ones(size(eTimes))*double(eThresh),'r*')
                    hold on
                    plot_if(bTimes*(alpha(it)),ones(size(eTimes))*double(eThresh)*1.2,'mx')
                    legend('beh Trigs');
                    plot(eTRun(eT0:end))
                end
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
        if strcmpi(figures,'plot')
            figure(sf(irec))
            title(sprintf('Shift estim + align, rec %s run %02d',info.rec(nrec).name, irun), 'FontSize', 10, 'FontWeight', 'bold')
        end
    end
end
save(fn.sess_info, 'info')
fprintf('trial shift info saved in %s\n',fn.sess_info)
    
end

function trials_match(mouse, sess, varargin)
%estimates the misalingment (in trials) of the h5 file of timestapms
%respect to the ephys recordings
%It was made for trigger signals in the beh files that count the time since
%the beginning of the trial.
%For use with the Trial Pin and back compatibility,adapted that the time tt
%does not add any other delay if 'runTrialStart' and 'TrPin' are the
%signals.
%global info


inPar=inputParser;
inPar.addRequired('mouse')
inPar.addRequired('sess')
inPar.addParamValue('figures','plot',@(x) ischar(x) && any(strcmpi(x,{'plot','noplot'})));
inPar.addOptional('rec','',@(x)iscell(x) || (ischar(x) && ~any(strcmpi(x,inPar.Parameters))) );

inPar.parse(mouse,sess, varargin{:});
figures=inPar.Results.figures;

fn = file_names(mouse, sess);
q = load(fn.sess_info);
info = q.info;
if ~isnumeric(sess)
    sess=str2num(sess);
end

fprintf('\n-----------------------------------------------------------\n');
fprintf('Estimating trial shift for mouse %s, session %s...\n',num2str(mouse),num2str(sess));

rec=inPar.Results.rec;
if nargin<3 || isempty(rec)
    fprintf('\t no record specified, going for all of them...\n');
    reclist = {q.info.rec.name};
else
    if iscell(rec)
        reclist =rec;
    else
        reclist={rec}
    end
end

eTrig_chan = 'TrPin';
bTrig_chan = 'runTrialStart';
for irec=1:numel(reclist)
    sf(irec)=figure;
    set(sf(irec),'NextPlot','add');
end

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

    
    for irun = 1:nrun
        run=info.rec(nrec).run(irun);
        runFig(irun)=figure();
        % e-phys trigger
        run_start = info.rec(nrec).run(irun).start;
        run_dur   = info.rec(nrec).run(irun).duration;
        
        %GET THE TRIGGERS
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
        
        %                     figure
        %                     plot(eTRun)
        %                     plot(diff(eTRun))
        %                     hold
        %                     plot(eSyncRun,ones(size(eSyncRun))*15000,'ro')
        %These are the trials for this run
        in=find( ([trial.run] == irun));
        bSyncRun = [bSync(in);bSync(in)+[trial(in).runTrialDur]];
        %%%% Now do the matching for the trials in this run
        %Have eSyncRun, bSync (all the trial starts for the rec), in (indexes of all the trials for this run);
        if irun==4
            disp('run 4')
        end
        [inMatched,eMatched]=match_trials_run(eSyncRun,bSyncRun,in)
        
        bMatched=bSyncRun(:,inMatched-in(1)+1);
        nn=length(inMatched);
        % estimation synchronization error
        error = nan;
        alpha = nan;
        
        %  and alinging triggers by the first trial
        % relative to beggining of run
        eT0=eMatched(1,1); % the t0 of the aligned ephys, segment, relative to beginning of the run
        bT0=bMatched(1,1);
        eTimes = eMatched(:,1:nn) - eMatched(1,1);
        bTimes = bMatched(:,1:nn) - bMatched(1,1);
        ind    = 1:nn;
        alpha  = (eTimes(1,:)/bTimes(1,:));
        error  = sqrt(norm(eTimes-alpha*bTimes))
        
        if strcmpi(figures,'plot')
            figure(sf(irec))
            subplot(nrun,3,3*irun-2)
            plot(eTimes,alpha*bTimes,'x')
            xlabel('ephys trig');
            ylabel('behav trig');
            title('trial_shift')
            subplot(nrun,3,3*irun-1)
            errorVec=sum(eTimes-alpha*bTimes,1);
            plot(errorVec, 'x')
            xlabel('trial');
            ylabel('error');
            figure(runFig(irun))
            clf
            plot(eTimes,ones(size(eTimes))*double(eThresh),'r*')
            hold on
            plot_if(bTimes*(alpha),ones(size(eTimes))*double(eThresh)*1.2,'mx')
            legend('beh Trigs');
            plot(eTRun(eT0:end))
        end
        title(sprintf('Shift estim + align, rec %s run %02d',info.rec(nrec).name, irun), 'FontSize', 10, 'FontWeight', 'bold')
        trAlign.eTrChan=eTrig_chan;
        trAlign.bTrField=bTrig_chan;
        trAlign.eT0Run = eT0; %start of the aligned ephys rel. to the run.
        trAlign.eTrigsT0 = eTimes; %ephys trig times rel to eT0
        trAlign.iTrialsRun = inMatched; %beh trials that have a matching ephys trial.
        
        %Align the trials for this run:
        %Get them zero start time by default
        for it=1:numel(inMatched)
            %there are as many beh trial (numel(inMatched)) as ephys
            %triggers (eMatched)
            % shift the beginning of the trial to the corresponding
            % ephys trigger
            %(it)th trigger start relative to the beginning of th rec
            t0=eMatched(1,it)+run_start;
            trial(inMatched(it)).start=t0;
            fprintf('%d    %d \n', inMatched(it), trial(inMatched(it)).start)
        end
    end  % for irun = 1:nrun:
    drawnow
    
    trAlign.frequencyCorrection = alpha;
    info.rec(nrec).run(irun).trAlign=trAlign;
    disp(info.rec(nrec).run(irun).trAlign)
    if strcmpi(figures,'plot')
        figure(sf(irec))
        title(sprintf('Shift estim + align, rec %s run %02d',info.rec(nrec).name, irun), 'FontSize', 10, 'FontWeight', 'bold')
    end
    
    save(fn.sess_info, 'info')
    fprintf('trial shift info saved in %s\n',fn.sess_info)
    save(fn.trial, 'trial')
    fprintf('trial structures saved in %s\n',fn.trial)
end




%_____________________________________________________%

    function [inMatched,eMatched] = match_trials_run(eT,bT,in)
        %get the segments within trials with zero start time
        %these might be trials not recorded in the behavior file or trials that didn't happen.
        %For each segment within two bad segments, find the shift and
        %match.
        %the segments:
        %assume all the trials from the beginning are good, find irruptions
        %of bad segments from there.
        %there will be transitions from good trials to bad trials, those
        %mark the end of "good" segments.
        %the transitions from bad to good marke the beginning of the next
        %"good" segments.
        
        trIn=trial(in); %all the trials from the beginning of the run

        all=zeros(1,numel(trIn));
        good=find(([trIn.runTrialDur]>0) & ([trIn.runTrialStart]>0));
        all(good)=1;
        goodToBad=find(diff(all)==-1); %last good trial of a chunk of good trials
        badToGood=find(diff(all)==1);  %last bad trial of a chunk of bad trials
        %if it starts bad:
        firstGood = 1;
        if ~(isempty(goodToBad) && isempty(badToGood))
        if (isempty(goodToBad) && ~isempty(badToGood)) || (~(isempty(goodToBad) || isempty(badToGood)) && badToGood(1)<goodToBad(1))
            firstGood = badToGood(1) +1;
            badToGood(1) = [];
        end
        end
        if numel(badToGood)==numel(goodToBad)
            goodSegments=[goodToBad;badToGood+1];
            goodSegments=reshape([firstGood goodSegments(:)' range(in)+1],2,[]);
        else
            goodSegments=[[firstGood badToGood];[goodToBad]];
        end

        
        inMatched=[];
        eMatched =[0;0];
        
        for i=1:numel(goodSegments(1,:))
            figHandles(i).trigs=figure(i+100+irun*10);
            clf;
            figHandles(i).alignment=figure(i+100+irun*10+numel(goodSegments(1,:)));
            clf;
            inS = goodSegments(1,i):goodSegments(2,i);
            % align the beh trials in the block (indexes inS)
            % wiht the segment of the eSR triggers that starts at the end
            % of the ones that have already been matched.
            % find the last matched
            firstUnmatched = find(eT(1,:)>eMatched(1,end),1);
            [eM,inM]=match_trials(eT(:, firstUnmatched: end),inS,trIn,figHandles(i));
            
            if strcmpi(figures,'plot')
                figure(figHandles(i).alignment)
                suptitle(sprintf('Shift estim + align, rec %s run %02d chunk %2d/%2d',info.rec(nrec).name, irun,i, numel(goodSegments(1,:))))
                figure(figHandles(i).trigs)
                title(sprintf('Shift estim + align, rec %s run %02d chunk %2d/%2d',info.rec(nrec).name, irun,i,numel(goodSegments(1,:))), 'FontSize', 10, 'FontWeight', 'bold')
            end
            eMatched=[eMatched eM];
            inMatched=[inMatched inM];
        end
        eMatched(:,1)=[];
        %get the indices back to the rec indexes
        inMatched=inMatched+in(1)-1;
    end

    
end

function [eMatched, iMatched]=match_trials(eSyncChunk,in,trial,figHandles)
% estimation synchronization error for each trial shift value
%whether the timing is relative to start of run is irrelevant
bSyncRun =[trial(in).runTrialStart];

% it will only work if the eSyncChunk array has at least three triggers.
if numel(eSyncChunk(1,:))<3 || numel(in)<3
    warning('Skipping chunk, it has less than 3 good trials');
    eMatched=[];
    iMatched=[];
    return
end

%jerry-rig for the subplotting

figures='plot';
nrun=1;
irun=1;
if nargin>3
    runFig=figHandles.trigs;
    sf=figHandles.alignment;
else
    runFig=figure;
    sf    =figure;
end
set(sf,'NextPlot','add');
set(runFig,'NextPlot','add');



trial_shift = [-5:5];
nts = numel(trial_shift);
error    = nan(1,nts);
alpha    = nan(1,nts);
err_min  = 1e10;
eThresh=1;


for it = 1:nts
    % defining the shifted trigger times
    %positive shift is when the first n trials of the beh were
    %missed
    %negative is when triggers dont have behavior and have to be skipped
    try
    if trial_shift(it) >= 0
        bTimes = bSyncRun(:,(1+trial_shift(it)):end);
        inShifted= in(1+trial_shift(it):end);
        eTimes = eSyncChunk;
    else
        bTimes = bSyncRun;
        eTimes = eSyncChunk(:,(1-trial_shift(it)):end);
        inShifted = in;
    end

    
    ne = length(eTimes);
    nb = length(bTimes);
    nn = min([ne, nb]);
    
    %  and alinging triggers by the first trial
    eT0=eTimes(1,1); % the t0 of the aligned ephys, segment, relative to beginning of the run
    bT0=bTimes(1,1);
    eTimes = eTimes(:,1:nn) - eT0;
    bTimes = bTimes(:,1:nn) - bT0;
    
    catch
        warning('Attempting an invalid trial shift (%2d)',trial_shift(it))
        continue
    end
    
    if (length(eTimes)<3)
        continue;
    end
    ind    = 1:nn;
    alpha(it) = (eTimes(1,:)/bTimes(1,:));
    error(it) = sqrt(norm(eTimes(1,:)-alpha(it)*bTimes(1,:)));
    
    if strcmpi(figures,'plot')
        figure(sf)
        subplot(nrun,3,3*irun)
        plot(trial_shift(it), error(it), 'o', 'Color', 'r', 'LineWidth', 2)
        xlabel('shift');
        ylabel('error sum');
        hold on
    end
    %
    if error(it)<err_min
        err_min = error(it);
        itmin   = it;
        if strcmpi(figures,'plot')
            figure(sf)
            subplot(nrun,3,3*irun-2)
            plot(eTimes,alpha(it)*bTimes,'x')
            xlabel('ephys trig');
            ylabel('behav trig');
            title(trial_shift(it))
            subplot(nrun,3,3*irun-1)
            errorVec=sum(eTimes(1,:)-alpha(it)*bTimes,1);
            plot(errorVec, 'x')
            xlabel('trial');
            ylabel('error');
            figure(runFig)
            clf;
            hold on;
            plot(eTimes,ones(size(eTimes))*double(eThresh),'*')
            hold on
            plot(bTimes*alpha(it),ones(size(bTimes))*double(eThresh)*1.2,'mx')
            ylim([0 1.5*double(eThresh)]);
            legend({'ephys Trigs','beh Trigs'});
            %plot(eTRun(eT0:end))
            %title(sprintf('Shift estim + align, rec %s run %02d',info.rec(nrec).name, irun), 'FontSize', 10, 'FontWeight', 'bold')
        end
        eMatched = eTimes+eT0; %ephys trig times rel to eT0
        iMatched = inShifted(ind); %indices of beh trials (trial) that have a matching ephys trial.
    end
    drawnow
end

end %function match_trials

function trial_alignment(mouse, sess, varargin)
%aligns the trials and the ephys file.

inPar=inputParser;
inPar.addRequired('mouse')
inPar.addRequired('sess')
inPar.addOptional('rec','',@(x)iscell(x) || (ischar(x) && ~any(strcmpi(x,inPar.Parameters))) );
inPar.parse(mouse,sess, varargin{:});

fn = file_names(mouse, sess);
q = load(fn.sess_info);
info = q.info;

    
fprintf('\n-----------------------------------------------------------\n');
fprintf('Aligning trials for mouse %s, session %s...\n',num2str(mouse),num2str(sess));
rec=inPar.Results.rec;
if nargin<3 || isempty(rec)
    fprintf('\t no record specified, going for all of them...\n');
    reclist = {q.info.rec.name};
else
    if iscell(rec)
        reclist =rec;
    else
        reclist={rec}
    end
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

function sniffs_to_trials(mouse, sess, varargin)
% the function reads the strucutre sniff.mat from xxx_sniff.mat file and
% sort the sniff parameters into individula trials, save the results tback
% to trial structure to the file xxx_trials.mat
% the mosue and sess must be specifiied. 
% rec is optional
% if rec is specified, it sort sniffs to trials for a specific recors,
% otehrwise, it reads the lsit of all records for a given session

inPar=inputParser;
inPar.addRequired('mouse')
inPar.addRequired('sess')
inPar.addOptional('rec','',@(x)iscell(x) || (ischar(x) && ~any(strcmpi(x,inPar.Parameters))) );
inPar.parse(mouse,sess, varargin{:});

fn = file_names(mouse, sess);
q = load(fn.sess_info);
info = q.info;

global trial

    fprintf('\n-----------------------------------------------------------\n');
    fprintf('Getting sniff times into trial structures for mouse %s, session %s...\n',num2str(mouse),num2str(sess));

    fn = file_names(mouse, sess);
    q = load(fn.sess_info);
    rec=inPar.Results.rec;
    if nargin<3 || isempty(rec)
        fprintf('\t no record specified, going for all of them...\n');
        reclist = {q.info.rec.name};
    else
        if iscell(rec)
            reclist =rec;
        else
            reclist={rec}
        end
    end

    for ir = 1:numel(reclist)
        fprintf('rec %s\n ================\n', reclist{ir})
        fn = file_names(mouse, sess, reclist{ir});
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
            fprintf('Rec %s has no sniff parabola fittings, leaving sniffParabolaZeroTimes field out\n',reclist{ir});
        end
        
        
        for it = 1:numel(trial)
            if isempty(trial(it).start)||(trial(it).start == 0)
                continue
            end
            
            t0 = trial(it).start;
            t1 = t0 - 5000;
            t2 = t0 + 6000;
            
            in_sn = (sn_t_zer(3,:) > t1)&(sn_t_zer(1,:) < t2);
            
            trial(it).sniffZeroTimes      = sn_t_zer([1,2], in_sn) - t0;
            if isfield(sniff,'t_zer_fit')
                trial(it).sniffParabZeroTimes = sn_t_zer_fit(:,in_sn)  - t0;
            end
            fprintf('%d    %d \n', it, sum(in_sn))
            
            % If odor trial, check if first inhalation comes before odor
            % stimulus shuts off. If not, mark trial with conc = -1.
            if diff([trial(it).odorTimes]) == 0
                continue
            else
                first_sn_in = find(trial(it).sniffZeroTimes(1,:) > trial(it).odorTimes(1),1);
                inh1_delay = trial(it).sniffZeroTimes(1,first_sn_in) - trial(it).odorTimes(1);
                bad_sniff_trials = inh1_delay > diff(trial(it).odorTimes);
                
                if isempty(inh1_delay) || bad_sniff_trials
                    sprintf('    bad trial: %i',it)
                    trial(it).odorConc = -1;
                end
            end
            
        end
        save(fn.trial, 'trial')
        fprintf('savedd sniffs in trial structure %s\n', fn.trial)
    end
end

function spikes_to_trials(mouse, sess, varargin)

inPar=inputParser;
inPar.addRequired('mouse')
inPar.addRequired('sess')
inPar.addParameter('newSpikes',false,@islogical);
inPar.addOptional('rec','',@(x)iscell(x) || (ischar(x) && ~any(strcmpi(x,inPar.Parameters))) );
inPar.parse(mouse,sess, varargin{:});

global trial

    fn = file_names(mouse, sess);
    q = load(fn.sess_info);
    
    fprintf('\n-----------------------------------------------------------\n');
    fprintf('Getting spike times into trial structures for mouse %s, session %s...\n',num2str(mouse),num2str(sess));
    rec=inPar.Results.rec;
    if nargin<3 || isempty(rec)
        fprintf('\t no record specified, going for all of them...\n');
        reclist = {q.info.rec.name};
    else
        if iscell(rec)
            reclist =rec;
        else
            reclist={rec}
        end
    end

    for ir = 1:numel(reclist)
        fprintf('rec %s\n ================\n', reclist{ir})
        if inPar.Results.newSpikes
            wrap_message(sprintf('Making spikes structures from sorting files'),'*');
            afta_da_sortin(mouse,sess,reclist{ir});
        end
        fn = file_names(mouse, sess, reclist{ir});
        q = load(fn.trial);
        trial = q.trial;
        
        q = load(fn.spikes);
        
        for iu = 1:numel(q.unit)
            sp{iu} = q.unit(iu).times;
        end
        
        
        for it = 1:numel(trial)
            trial(it).spikeTimes=[];
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
