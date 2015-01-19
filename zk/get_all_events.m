function [allEvents, aev_m] = get_all_events(mouse, sess, rec,irun,varargin)
%get all the events on a run
%events are extracted from the corresponding channel of a binary file
%it returns a matlab structure.

%trialNumbers(:,i)=[onTime
%                   words
%                   ]
%where words is a structure containing the bytes.
event = struct('event',[],'type',[],'chanId',[],'evtFcn',[]);
inputPar = struct('par',[],'default',[],'validation',[]);

eventsList = [];

%add all the events
%trial number
event.event  = 'trialNumber';
event.type   = 'serial';
event.chanId = 'trNum';
event.evtFcn = @get_trial_numbers;
eventsList   = [eventsList event];

% trial pin
event.event  = 'trialPin';
event.type   = 'digital';
event.chanId = 'trPin';
event.evtFcn = @get_analog_events;
eventsList   = [eventsList event];

% final valve pin
event.event  = 'finalValve';
event.type   = 'digital';
event.chanId = 'FVPin';
event.evtFcn = @get_analog_events;
eventsList   = [eventsList event];

% laser
event.event  = 'laser';
event.type   = 'semiDigital';
event.chanId = 'Laser';
event.evtFcn = @get_analog_events;
eventsList   = [eventsList event];



inPar=inputParser;
inPar.addRequired('mouse')
inPar.addRequired('sess')
inPar.addRequired('rec',@ischar)
inPar.addRequired('irun',@isscalar)
inPar.addParamValue('figures','plot',@(x) ischar(x) && any(strcmpi(x,{'plot','noplot'})));
inPar.addParamValue('aqSystem','JFRC',@(x) ischar(x) && any(strcmpi(x,{'JFRC','intan'})));
inPar.addParamValue('eventsList',eventsList,@(x) isstruct(x) && all(strcmp(fields(x),fields(event))));
inPar.addParamValue('getSniff',false,@islogical);
inPar.addParamValue('getSpikes',false,@islogical);
% serial data parameters
inPar.addParamValue('serial',{300,8,0,1,1000},@(x) iscell(x) && numel(x)==5);
inPar.addParamValue('wordLength',2,@(x) isscalar(x));
inPar.addParamValue('samplingFrequency',19385,@(x) isscalar(x))
inPar.addParamValue('threshold',0.1,@(x) isscalar(x));
inPar.parse(mouse,sess, rec,irun,varargin{:});
figsOn=inPar.Results.figures;

if ~isnumeric(sess)
    sess=str2num(sess);
end

fprintf('\nGetting all events for %s, run %02d...\n',[mouse '_' num2str(sess,'%02d') 'rec'],irun);
fprintf('--------------------------------------------------\n')
%get the info for the rec and for the run
fn = file_names(mouse, sess,rec);
q = load(fn.sess_info);
recInfo = q.info.rec(strcmpi(rec,{q.info.rec.name}));
runInfo = recInfo.run([recInfo.run.num]==irun);




%get the channel from the JFRC system:
%get the file:
runRawFile= fullfile(fn.fold_rd_sess,runInfo.ephys_data);
h5file    = fullfile(fn.fold_pr_sess,[fn.basename_an 'evt.h5']);
%read from the info the filename from the run and get the binary channel
%create the structure with all the events for this mouse, sess, rec, run
allEvents   = struct();
foundEvents = ones(1,numel(eventsList));
for iEvent = 1:numel(eventsList)
    %get the info on the channel
    fprintf('\n** Event %d/%d (%s, type %s) **\n',iEvent,numel(eventsList),eventsList(iEvent).event,eventsList(iEvent).type);
    evt    = eventsList(iEvent);
    eventChanNum = find(strcmpi(evt.chanId,{recInfo.chan.name}));
    if ~isempty(eventChanNum)
        chInfo = recInfo.chan(eventChanNum);
        %read the raw data
        %rawStream = [];
        rawStream = read_raw_stream(chInfo);
        %make the stream into events
        allEvents.(evt.event) = evt.evtFcn(rawStream);
    else
        warning('Did not find event %s (looking for channel name %s)',evt.event,evt.chanId);
        foundEvents(iEvent) = 0;
    end
end
eventsList=eventsList(find(foundEvents));
fprintf('\nDone getting events (%d event streams found)\n',numel(eventsList));
fprintf('--------------------------------------------------\n')
% save the structure with all the events
%write into the events file
%write_h5(allEvents.(evt.event),evt,h5file);
aev_m = write_m(allEvents,eventsList,fn.evt_m,runInfo);
%if the option was set, get other events (sniff, spikes)
if inPar.Results.getSniff
    sniffs2events(mouse,sess,rec,irun);
end
if inPar.Results.getSpikes
    spikes2events(mouse,sess,rec,irun);
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function stream = read_raw_stream(eventChan)
        %read a channel with recorded events
        %eventChan = structure with channel info
        %stream    = structure with fields 
        %   data = the raw data in the channel
        %   int2Volt = conversion factor from integer value to voltage
        %   threshold = threshold (for on/off extraction) in integers
        if strcmpi(inPar.Results.aqSystem,'JFRC')
            fprintf('Reading events channel (%s)...',eventChan.name);
            runFID=fopen(runRawFile,'r');
            stream.data = read_analog_channel(runFID,eventChan.num_rd,recInfo.nChan_rd);
            fclose(runFID);
            fprintf('done.\n');
            
            %it is easier to work with integers, so get the Voltage values to integers
            chanGain = 1;
            if isfield(eventChan,'gain') && ~isempty(eventChan.gain)
                chanGain = eventChan.gain;
            end
            volt2int=32768/(recInfo.rangeMax); % in V
            stream.int2Volt = 1./volt2int;
            stream.chanGain = chanGain;
            
            %the minimum detectable peak is 0.5v
            stream.threshold=ceil(inPar.Results.threshold*volt2int*chanGain);
        else
            error('Dont know how to handle aqcuisition system %s yet',inPar.Results.aqSystem);
        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~
    function analogEvents = get_analog_events(stream)
        %get the analog events
        %input is a stream structure with fields:
        % data      : the stream
        % int2Volt  : factor from int 2 volt value
        % threshold : threshold for onset of event
        % retunrs analogEvents structure array with fields:
        % on
        % off
        % value (int)
        fprintf('Extracting analog events from channel...')
        eventOnOff=zeros(size(stream.data));
        eventOnOff(stream.data>stream.threshold)=1;
        
        eventOnSample  = find(diff(eventOnOff)==1);
        eventOffSample = find(diff(eventOnOff)==-1);
        
        if eventOffSample(1)<eventOnSample(1)
            eventOffSample(1)=[];
        end
        if eventOnSample(end)>eventOffSample(end)
            eventOnSample(end)=[];
        end
        
        if ~(numel(eventOnSample)==numel(eventOffSample))
            error('Inconsistent number of on/off events: check threshold, maybe?');
        end
        
        %get the mean values of every ON segment and make a struct array of
        %event = struct('on',[],'off',[],'value,[])
        analogEvents = arrayfun(@(x,y) struct('on',x,'off',y, 'value',mean(stream.data(x+1:y-1))*stream.int2Volt/stream.chanGain),eventOnSample,eventOffSample);
        fprintf(' done.\n');
    end
    % ~~~~~~~~~~~~~~~~~~~~~~
    function trialNumbers = get_trial_numbers(stream) 
        %input is a stream structure with fields:
        % data      : the stream
        % int2Volt  : factor from int 2 volt value
        % threshold : threshold for onset of event
        % retunrs analogEvents structure array with fields:
        % on
        % off
        % value (int)
        %make a serialdata object, fill it with the stream, and get the
        %numbers
        %make the serial object and fill it
        fprintf('Extracting serial numbers from channel...')
        serialPars = {'serial','wordLength','samplingFrequency'};
        varSerial = cell(1,numel(serialPars)*2);
        for ip=1:numel(serialPars)
            varSerial{ip*2-1}   = serialPars{ip};
            varSerial{ip*2} = inPar.Results.(serialPars{ip});
        end
        serialNumbers = serialdata(varSerial{:});
        serialNumbers.put_stream(stream.data);
        %get the trial numbers
        [stamps, numbers] = serialNumbers.get_numbers();
        trialNumbers = arrayfun(@(x,y) struct('on',x,'value',y),stamps,numbers);
        fprintf(' done.\n');
    end
    % ~~~~~~~~~~~~~~~~~~~~~~
    function events = write_m(allEvents,eventsMeta,evFileName,runInfo)
        % open the event file for the record (%create it if it doesnt
        % exist),
        % save the events and meta under an field 'xx' (two digits for the number of run)
        [evDir,evFile,evExt] = fileparts(evFileName);
        runFieldName     = sprintf('run_%02d',runInfo.num);
        fprintf('\n - Saving events structure of run %s in matlab file %s%s...\n\t ',runFieldName,evFile,evExt);
        % if the file exists, open it, if not, make it
        evFound = dir(evFileName);
        if ~isempty(evFound)
            fprintf('Found an events matlab file for this rec. I will open it and re-write the events for this run.\n');
            events = load(evFileName);
            if isfield(events,runFieldName)
                events=rmfield(events,runFieldName);
            end
        else
            fprintf('No previous matlab event file for this rec. Making a new one.\n');
            events = struct();
        end
        
        %add the event structure for this run.
        for iev=1:numel(eventsMeta)
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
        %save
        if ~exist(evDir,'dir')
            mkdir(evDir)
        end
        save(evFileName,'-struct','events');
        fprintf(' - done writing events file in pr_data folder.\n');
        
    end

    % ~~~~~~~~~~~~~~~~~~~~~~
    function write_h5_file(horizEvents,fileName,runInfo)
        %write events for the run to the h5 file
        %checks if file exists. if it doesnt, it creates it.
        %opens the file and replaces the group for the run (deletes if
        %exist)
        %creates the group for the run and makes a table (vector) with
        %everyone of the events in the horizontalized events structure.
        if ~exist(fileName,'file')
            fprintf('No previous hdf5 event file for this rec. Making a new one.\n');
            fcpl = H5P.create('H5P_FILE_CREATE');
            fapl = H5P.create('H5P_FILE_ACCESS');
            fid  = H5F.create(fileName,'H5F_ACC_TRUNC',fcpl,fapl);
        else
            fprintf('Found an events hdf5 file for this rec. I will open it and re-write the events for this run.\n')
        end
        %create the group for the run if it doesnt exist
        h5i=h5info(fileName);
        runGroupName=['run_' num2str(runInfo.num,'%02d')];
        if any(strcmpi(runGroupName,{h5i.Group
        
        
        %close the file
        H5F.close(fid);
    end

    % ~~~~~~~~~~~~~~~~~~~~~~
    function Y = read_analog_channel(fid, ich, nch)
        fseek(fid, 2*(ich-1), 'bof');
        Y = fread(fid, inf, 'int16', (nch-1)*2);
    end
    % ~~~~~~~~~~~~~~~~~~~~~~
end