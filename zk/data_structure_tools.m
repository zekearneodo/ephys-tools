%Set of tools for extracting events and making trial structures.
%
%
function ds = data_structure_tools()
global ds;
[~,computerName]=system('hostname');
if ~strcmp(strtrim(computerName),'flipper')
    includePath=fullfile(fileparts(pwd),'current','include');
    addpath(includePath);
end
    ds.events_lookup       = @events_lookup; %lookup a list of events in the corresponding trial in the voyeur table
    ds.match_trial_numbers = @match_trial_numbers; %get trial numbers and match them with trial pin events
    ds.make_fv_table       = @make_fv_table; % make table for final valve events
    ds.make_laser_table    = @make_laser_table; % make table for final valve events
    ds.make_tables         = @make_tables;
    ds.get_analog_events   = @get_analog_events; %get occurrences of an analog event 
    ds.get_trial_numbers   = @get_trial_numbers; % get the trial numbers
    ds.get_data_stream     = @get_data_stream; %get stream of data from binary file
    ds.get_voyeur_table    = @get_voyeur_table; %get a voyeur table (or fields of it)
    ds.read_voyeur_fields  = @read_voyeur_fields; %read fields for an in voyeur file
    ds.voyeur_definitions  = @voyeur_definitions; %get a structure with parameter name translations for an event 
    ds.get_info            = @get_info; %gets info structures for rec, run
    ds.check_list          = @check_list; %checks a list for repeated or skipped positions
    
    ds.h52m_table          = @h52m_table; %reformats an h5 table into matlab data types
end

function evt = make_tables(mouse,sess,rec,run)
%make all the event tables for a run
%these tables will be then gathered and put together in a table for each
%rec

%where words is a structure containing the bytes.
event = struct('event',[],'type',[],'chanId',[],'evtFcn',[]);
inputPar = struct('par',[],'default',[],'validation',[]);

%%% DEFINITIONS OF EVENTS AND PARAMETERS
% this will eventually go in a file for the session or the experiment
eventsList = [];

% trial pin
event.name     = 'trialPin';
event.type     = 'digital';
event.chanId   = 'trPin';
event.evtFcn   = @get_analog_events;
event.tableFcn = @check_trial_pin;
%eventsList     = [eventsList event];

% final valve pin
event.name     = 'finalValve_1';
event.type     = 'digital';
event.chanId   = 'FVPin';
event.evtFcn   = @get_analog_events;
event.tableFcn = @make_fv_table;
eventsList     = [eventsList event];

% laser
event.name     = 'laser_1';
event.type     = 'semiDigital';
event.chanId   = 'Laser';
event.evtFcn   = @get_analog_events;
event.tableFcn = @make_laser_table;
eventsList     = [eventsList event];


%get and match the trial numbers with the trial pins
trialNumbers = match_trial_numbers(mouse,sess,rec,run);

% for every event, look it up within its corresponding trial,
% and make the table with all its properties
for ie = 1:numel(eventsList);
    event = eventsList(ie);
    eventsTrials = events_lookup(event,trialNumbers,mouse,sess,rec,run);
end


end

function tn = match_trial_numbers(mouse,sess,rec,run)
%gets the time stamps of trial pins (on/off)
%gets the trial number that corresponds to that pin
%makes a table of the trial number and its correspondign pin
%get the trial numbers
tn.name   = 'trialNumber';
tn.type   = 'serial';
tn.chanId = 'trNum';

tn.events = get_trial_numbers(mouse,sess,rec,run,'figures','noplot');

%get the trial pins
tp.name   = 'trialPin';
tp.type   = 'digital';
tp.chanId = 'trPin';
tp.events = get_analog_events(tp.chanId,mouse,sess,rec,run,'figures','noplot');

%assign to every element in trial_number the immediate next trialPin
% 2 steps:
% pair every trNumber with the first next trPin
% look for n-uplicates and discard the first
%get the first pin that comes after a trial number.
%there can be a degenaracy: several trial numbers have the same pin.
%this means a trial number was sent in between trial pins, in which case
%we want to keep the closest (the latest) trial number to that pin
%the pins after every trial number
numStartPins = arrayfun(@(x) find(tp.events(1,:)>x,1),[tn.events.on]);
tStartPins   = tp.events(1,numStartPins);
tEndPins     = tp.events(2,numStartPins);
%look for errors in the list
er = check_list(numStartPins);
if ~isempty(find(er==1, 1))
    warning('Some trial pins did not have a matching trial number');
end
if ~isempty(find(er==2, 1))
    warning('Some trial pins have several trial numbers in between');
    %find the repeated numbers and un-match the further trial numbers from
    %a degenerate trial pin
    repTP=unique(numStartPins(diff(numStartPins)==0));
    fprintf('Solving %d degenerate trial pins:\n',numel(repTP));
    for pin=repTP
         %gather the indices of the trial numbers that share the same pin
         %get the degenerate pin event numbers
         degTId=find(numStartPins==pin);
         rightTN = degTId(find([tn.events(degTId).on]<tp.events(1,pin),1,'last'));
         badTId  = degTId(~(degTId==rightTN));
         numStartPins(badTId) = nan;
         tStartPins(badTId)   = nan;
         tEndPins(badTId)     = nan;
    end
end

%todo: check trial duration consistence

tn.table.trialNumber = [tn.events.value];
tn.table.trialOn     = [tn.events.on];
tn.table.trialPin    = numStartPins;
tn.table.trialStart  = tStartPins;
tn.table.trialEnd    = tEndPins;
%tn.meta = %some metadata that will be added when save in h5 format

end

function get_laser(mouse,sess,rec,irun)
%gets laser time stamps and v value
%gets the trial number that corresponds to that laser
%gets the laser parameter that correspond to that laser event


end

function [allEvents, aev_m] = make_events_tables(mouse, sess, rec,irun,varargin)
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
inPar.addParamValue('samplingFrequency',19531.25,@(x) isscalar(x))
inPar.addParamValue('threshold',0.25,@(x) isscalar(x));
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

    end

    % ~~~~~~~~~~~~~~~~~~~~~~
    function Y = read_analog_channel(fid, ich, nch)
        fseek(fid, 2*(ich-1), 'bof');
        Y = fread(fid, inf, 'int16', (nch-1)*2);
    end
    % ~~~~~~~~~~~~~~~~~~~~~~
end

function [onEvents]=get_analog_events(chanName,mouse, sess, rec,irun,varargin)
%turns an analog channel with events (Hi/Lo) into a semi-digital:
%onEvents(:,i)=[onTime
%               offTime
%               onValue]
%where on/off Times are in sample number and onValue is the mean value in
%Volts.
%through the on segment.
%onValue is computed using the range and the number of bits of the
%acquisition system (taken from the info file) and the gain (individual for
%each channel, default is 1.

inPar=inputParser;
inPar.addRequired('chanName',@ischar);
inPar.addRequired('mouse')
inPar.addRequired('sess')
inPar.addRequired('rec',@ischar)
inPar.addRequired('irun',@isscalar)
inPar.addParameter('figures','plot',@(x) ischar(x) && any(strcmpi(x,{'plot','noplot'})));
inPar.addParamValue('aqSystem','JFRC',@(x) ischar(x) && any(strcmpi(x,{'JFRC','intan'})));
inPar.addParamValue('threshold',0.25,@(x) isscalar(x)); %threshold in volts
inPar.addParamValue('chanGain',1,@(x) isscalar(x));

inPar.parse(chanName,mouse,sess, rec,irun,varargin{:});
figsOn=strcmpi('plot',inPar.Results.figures);
eventChanName=inPar.Results.chanName;

%get the info for the rec and for the run
fn = file_names(mouse, sess,rec);
q = load(fn.ss_sess_info);
recInfo = q.info.rec(strcmpi(rec,{q.info.rec.name}));
runInfo = recInfo.run([recInfo.run.num]==irun);

if ~isnumeric(sess)
    sess=str2num(sess);
end


%get the channel from the JFRC system:
%get the file:
runRawFile=fullfile(fn.fold_rd_sess,runInfo.ephys_data);

%read from the info the filename from the run and get the binary channel
if strcmpi(inPar.Results.aqSystem,'JFRC')
    eventChan=recInfo.chan(strcmpi(eventChanName,{recInfo.chan.name}));
    fprintf('Reading analog events channel (%s)...',eventChan.name);
    runFID=fopen(runRawFile,'r');
    eventData=read_analog_channel(runFID,eventChan.num_rd,recInfo.nChan_rd);
    fclose(runFID);
    fprintf('done.\n');
    
    %it is easier to work with integers, so get the Voltage values to integers
    volt2int=32768/(recInfo.rangeMax*inPar.Results.chanGain); % in V
    int2volt=1./volt2int;
    
    %the minimum detectable peak is 0.5v
    threshold=ceil(inPar.Results.threshold*volt2int);
end

eventOnOff=zeros(size(eventData));
eventOnOff(eventData>threshold)=1;

eventOnSample  = find(diff(eventOnOff)==1);
eventOffSample = find(diff(eventOnOff)==-1);

if eventOffSample(1)<eventOnSample(1)
    eventOffSample(1)=[];
end
if eventOnSample(end)>eventOffSample(end)
    eventOnSample(end)=[];
end

if numel(eventOnSample)==numel(eventOffSample)
    eventSample=[eventOnSample';eventOffSample'];
else
    error('Inconsistent number of on/off events: check threshold, maybe?');
end

%get the mean values of every ON segment
eventAmp   =arrayfun(@(x,y) mean(eventData(x+1:y-1)),eventSample(1,:),eventSample(2,:));

if figsOn
    figure
    hold on
    plot(eventData)
    %plot(laserOnOff*40000,'m')
    plot(eventOnSample,eventAmp,'r*')
    plot(eventOffSample,eventAmp,'go')
end

onEvents=[eventSample;eventAmp*int2volt];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function Y = read_analog_channel(fid, ich, nch)
        fseek(fid, 2*(ich-1), 'bof');
        Y = fread(fid, inf, 'int16', (nch-1)*2);
    end

end

function table = events_lookup(ev,tn,mouse,sess,rec,irun)
%easier way to lookup trial numbers for an event
%  Go through the trial numbers that have a good trial pin
%  get the closest event after that pin (whithin the pin on and off)
%  lookup that pin in the trialnumber table and get the trial number
% ev.name   : name of the event
% ev.type   : type of event (for extraction of the pins)
% ev.chanId : id of the channel where the event is
% ev.evtFcn : pointer to the function used to extract timestamps and values
% Example:
% event.name   = 'finalValve_1';
% event.type   = 'digital';
% event.chanId = 'FVPin';
% event.evtFcn = @get_analog_events;
% get the events

events = get_analog_events(ev.chanId,mouse,sess,rec,irun,'figures','noplot');
% go through all the trial numbers that have a good trial pin and find the
% events whithin that trial
eventsTrials = [];
for itp = find(~isnan(tn.table.trialPin))
    tPinStamps = [tn.table.trialStart(itp);tn.table.trialEnd(itp)];
    trialNumber = tn.table.trialNumber(itp);
    if any(isnan(tPinStamps)) || ~(diff(tPinStamps)>0)
        warning('Wrong timestamps in trial %d\n',trialNumber);
        continue
    end
    % lookup the event(s) within the trial
    % (get the rows in the event table, then check them, then append them
    % to the table)
    evIdx = find((events(1,:)>=tPinStamps(1)) & (events(2,:) <= tPinStamps(2)));
    if ~isempty(evIdx)
        eventsTrials = [eventsTrials; {trialNumber, evIdx}];
    end
end
% this yelds a list of {trial number, index of event pins within % trial}
% now send the events and the eventsTrials to the tableFcn that corresponds
% to that event and make the table
table = ev.tableFcn(ev,events, eventsTrials,mouse,sess,rec,irun);
% now do some table to h5 something

end

function to=make_fv_table(ev,events,evTrials,mouse,sess,rec,irun)
 %receives eventsTrials (cell array of {trialNumber, eventIndex}
 %with trial Numbers matching timestamps.
 %make a table with all the data on every event.
 %the evTrials are checked to match an event.
 %just fill the table with the data on that fv open event.
 % go through the trials, check the trial has an open valve, and fill the
 % table.
 [rawTrials, rawTable] = read_voyeur_fields(ev.name, mouse,sess,rec,irun);
 %filter only the trials that come with clean corresponding events from
 %evTrials
 trials = [];
 
 for iEv = 1:length(evTrials)
     thisTrial    = rawTrials([rawTrials.trialNumber]==evTrials{iEv,1});
     thisEventIdx = evTrials{iEv,2};
     %check that there is only one event (one single final valve opening)
     if numel(thisEventIdx)>1
         warning('More than one opening final valve event in trial %d',thisTrial.trialNumber)
         continue
     end
     %if its ok, add the data of the event to the trial fields
     thisTrial.on  = events(1,thisEventIdx);
     thisTrial.off = events(2,thisEventIdx);
     trials = [trials thisTrial];
 end
 to.trials = trials;
end

function to=make_laser_table(ev,events,evTrials,mouse,sess,rec,irun)
 %receives eventsTrials (cell array of {trialNumber, eventIndex}
 %with trial Numbers matching timestamps.
 %make a table with all the data on every event.
 %the evTrials are checked to match an event.
 %just fill the table with the data on that fv open event.
 % go through the trials, check the trial has an open valve, and fill the
 % table.
 [rawTrials, rawTable] = read_voyeur_fields(ev.name, mouse,sess,rec,irun);
 %filter only the trials that come with clean corresponding events from
 %evTrials
 trials = [];
 
 for iEv = 1:length(evTrials)
     thisTrial    = rawTrials([rawTrials.trialNumber]==evTrials{iEv,1});
     thisEventIdx = evTrials{iEv,2};
     %check that there is only one event (one single laser pulse)
     if numel(thisEventIdx)>1
         warning('More than one laser event in trial %d',thisTrial.trialNumber)
         continue
     end
     %if its ok, add the data of the event to the trial fields
     thisTrial.on  = events(1,thisEventIdx);
     thisTrial.off = events(2,thisEventIdx);
     thisTrial.v   = events(2,thisEventIdx);
     trials = [trials thisTrial];
 end
 to.trials = trials;
end

function [trialNumbers]=get_trial_numbers(mouse, sess, rec,irun,varargin)
%reads a channel with the serial-transmitted trial numbers:
%trialNumbers(:,i)=[onTime
%                   words
%                   ]
%where words is a structure containing the bytes.
inPar=inputParser;
inPar.addRequired('mouse')
inPar.addRequired('sess')
inPar.addRequired('rec',@ischar)
inPar.addRequired('irun',@isscalar)
inPar.addParamValue('figures','plot',@(x) ischar(x) && any(strcmpi(x,{'plot','noplot'})));
inPar.addParamValue('aqSystem','JFRC',@(x) ischar(x) && any(strcmpi(x,{'JFRC','intan'})));
inPar.addParamValue('threshold',0.25,@(x) isscalar(x));
inPar.addParamValue('chanGain',1,@(x) isscalar(x));
inPar.addParamValue('chanName','trNum',@(x) ischar(x) );
% serial data parameters
inPar.addParamValue('serial',{300,8,0,1,1000},@(x) iscell(x) && numel(x)==5);
inPar.addParamValue('wordLength',2,@(x) isscalar(x));
inPar.addParamValue('samplingFrequency',19531.25,@(x) isscalar(x));

inPar.parse(mouse,sess, rec,irun,varargin{:});
figsOn=inPar.Results.figures;
eventChanName=inPar.Results.chanName;

%get the info for the rec and for the run
fn = file_names(mouse, sess,rec);
q = load(fn.ss_sess_info);
recInfo = q.info.rec(strcmpi(rec,{q.info.rec.name}));
runInfo = recInfo.run([recInfo.run.num]==irun);

if ~isnumeric(sess)
    sess=str2num(sess);
end

fprintf('Getting trial numbers for %s, run %02d...\n',fn.basename_an(1:end-1),irun);

%get the channel from the JFRC system:
%the channel:
trnChan=recInfo.chan(strcmpi(eventChanName,{recInfo.chan.name}));
%get the file:
runRawFile=fullfile(fn.fold_rd_sess,runInfo.ephys_data);

%make the serialdata object
serialNumbers = serialdata('serial',inPar.Results.serial,'wordLength',inPar.Results.wordLength,'samplingFrequency',inPar.Results.samplingFrequency);


%read from the info the filename from the run and get the binary channel
if strcmpi(inPar.Results.aqSystem,'JFRC')
    serialNumbers.read_stream_bin(runRawFile,trnChan.num_rd,recInfo.nChan_rd,'int16');
else
    error('Dont know how to handle %s acquisition system',inPar.Results.aqSystem);
end

%get the trial numbers
[stamps, numbers] = serialNumbers.get_numbers();
trialNumbers = arrayfun(@(x,y) struct('on',x,'value',y),stamps,numbers);

end

function [trials, table] = get_voyeur_table(mouse,sess,rec,run,varargin)
%reads the table in a voyeur-generated behavior file for a run.
% Assumes there is only one table and that table is the trial data
% Doesn't care about the streams
% returns:
% trials : array of trial structures with data in mlab format (double and char)
% table  : h5-formatted table, as it was read.
inPar=inputParser;
inPar.addRequired('mouse')
inPar.addRequired('sess')
inPar.addRequired('rec',@ischar)
inPar.addRequired('run',@isscalar)
inPar.addParameter('selectedFields',[], @iscell);

inPar.parse(mouse,sess, rec,run,varargin{:});


[~, runInfo] = get_info(mouse,sess,rec,run);
fn = file_names(mouse,sess,rec);

fnam    = fullfile(fn.fold_rd_sess, runInfo.behav_data);
hinfo   = h5info(fnam);

n_tr = numel(hinfo.Groups);
table_name = hinfo.Datasets.Name;
table = h5read(fnam,['/', table_name]);
ff = fieldnames(table);

if ~isempty(inPar.Results.selectedFields)
    keepFields = inPar.Results.selectedFields;
    %check that all the fields exist
    checked = cellfun(@(x) any(strcmp(x,ff)),keepFields);
    if any(~checked)
        warning('Some fields requested not found in the voyeur table');
    end
    %remove all the fields that do not appear in keepFields
    rmf = ~cellfun(@(x) sum(strcmp(x,keepFields)),ff);
    fieldsToRemove = ff(find(rmf));
    table = rmfield(table,fieldsToRemove);
    ff(find(rmf))=[];
end

%Make a structure array of this table
trials = struct();

for kf = 1:numel(ff)
    for it = 1:n_tr
        if ischar(table.(ff{kf}))
            trials(it).(ff{kf}) = deblank(table.(ff{kf})(:,it)');
        else
            trials(it).(ff{kf}) = double(table.(ff{kf})(it,1));
        end
    end
end


end

function stream = get_data_stream(chanName, mouse,sess,rec,irun,varargin)
%read a stream of data from a binary file
        %read a channel with recorded events
        %eventChan = structure with channel info
        %stream    = structure with fields 
        %   data = the raw data in the channel
        %   int2Volt = conversion factor from integer value to voltage
        %   threshold = threshold (for on/off extraction) in integers
        
        inPar=inputParser;
        inPar.addRequired('chanName',@ischar);
        inPar.addRequired('mouse')
        inPar.addRequired('sess')
        inPar.addRequired('rec',@ischar)
        inPar.addRequired('irun',@isscalar)
        inPar.addParameter('threshold',0.25,@(x) isscalar(x));
        inPar.addParameter('chanGain',1,@(x) isscalar(x));
        
        inPar.parse(chanName,mouse,sess, rec,irun,varargin{:});
        
        [recInfo,runInfo] = get_info(mouse,sess,rec,irun);
        streamChan=recInfo.chan(strcmpi(chanName,{recInfo.chan.name}));
        if isempty(streamChan)
            error('Channel %s not found in the dataset',chanName);
        end
        fprintf('Reading events channel (%s)...',streamChan.name);
        fn = file_names(mouse, sess,rec);
        runRawFile=fullfile(fn.fold_rd_sess,runInfo.ephys_data);
        runFID=fopen(runRawFile,'r');
        stream.data = read_analog_channel(runFID,streamChan.num_rd,recInfo.nChan_rd);
        fclose(runFID);
        fprintf('done.\n');
      
        %it is easier to work with integers, so get the Voltage values to integers
        chanGain = 1;
        if isfield(streamChan,'gain') && ~isempty(streamChan.gain)
            chanGain = streamChan.gain;
        end
        volt2int=32768/(recInfo.rangeMax); % in V
        stream.int2Volt = 1./volt2int;
        stream.chanGain = chanGain;
        
        %the minimum detectable peak is 0.5v
        stream.threshold=ceil(inPar.Results.threshold*volt2int*chanGain);
        
    % ~~~~~~~~~~~~~~~~~~~~~~
    function Y = read_analog_channel(fid, ich, nch)
        fseek(fid, 2*(ich-1), 'bof');
        Y = fread(fid, inf, 'int16', (nch-1)*2);
    end
    % ~~~~~~~~~~~~~~~~~~~~~~

end

function [recInfo, runInfo] = get_info(mouse,sess,rec,irun)
%get the info for the rec and for the run
fn = file_names(mouse, sess,rec);
q = load(fn.sess_info);
recInfo = q.info.rec(strcmpi(rec,{q.info.rec.name}));
runInfo = recInfo.run([recInfo.run.num]==irun);
end

function errorCode = check_list(list)
% checks on a list of numbers
% returns a list of codes:
% [] if everytihg is ok
% 1: skipped numbers
% 2: repeated numbers
errorCode = [];

%check if there were skipped numbers
if sum(diff(list)>1)
    errorCode = [errorCode 1];
    warning('Skipped numbers in the list');
end

if sum(diff(list)<1)
    errorCode = [errorCode 2];
    warning('Repeated numbers in the list');
end
    
end

function [trials, table] = read_voyeur_fields(evName,mouse,sess,rec,run)
% reads a voyeur table and translates the fields for the right type of
% event
%get the list of parameters saved and how their names translate to voyeur
%table
vd=voyeur_definitions(evName,mouse,sess,rec,run);
%vd.trialFields %for the trial
%vd.evtFields %for the event

%get all those fields from the table and convert the table to matlab
%datatypes
fieldsList = [ vd.trialFields ; vd.evtFields];
[trials, table]= get_voyeur_table(mouse,sess,rec,run,'selectedFields',fieldsList(:,2));
table=h52m_table(table);

%translate the names of the fields into the program names
for jf=1:length(fieldsList)
    oldField = fieldsList{jf,2};
    newField = fieldsList{jf,1};
    if strcmp(oldField,newField)
        continue
    else
        [trials.(newField)] = trials.(oldField);
        [table.(newField)]  = table.(oldField);
        trials = rmfield(trials, oldField);
        table  = rmfield(table, oldField);
    end
end

end

function vd = voyeur_definitions(evName,mouse,sess,rec,run)
%it uses default voyeur_definitions for name translations, unless there is
%a voyeur_definitions function in the mouse folder or addressed in the 
%here it checks if its indicated that it has to use other definitions
%vd_check
%(if it finds them; loads them and returns. if it doesnt, goes ahead)
%definitions:
% new columns:
%it returns vd structure with
%  trialFields : fields describing the trial 
%  evtFields   : fields describing the event named evName
% each variable xFields is a cell array of {'nameInDataStructure',
% 'nameInVoyeur'};

%description of trials
vd.trialFields = [ {'trialNumber','trialNumber'}
                {'trialComplete','trial_complete'}
                {'trialDuration','trialdur'}
                {'trialStart', 'starttrial'}
                {'trialEnd', 'endtrial'}
                {'trialIti', 'iti'}
                {'trialType','trialtype'}
                {'stimType','stimtype'}
                {'stimDescription','stim_desc'}
                {'noSniff','no_sniff'}
                {'paramsGot','paramsgottime'}
                {'result','result'}
                {'sniffDelay','sniffmaxdelay'}
                {'inhOnset','inh_onset'}
                ];
            
if strfind(evName,'finalValve')
    vd.evtFields = extract_finalValve_pars(evName);
    return
end

if strfind(evName,'laser')
    vd.evtFields = extract_laser_pars(evName);
    return
end
   

%~~~~ the functions for getting the tables

    function evtFields=extract_finalValve_pars(evName)
        evNameSplit=strsplit(evName,'_');
        evNum = (evNameSplit{end});
        %generate the fields to get from the voyeur_table
        evtFields = [  {['AirFlow_' evNum],['AirFlow_' evNum]}
                        {['NitrogenFlow_' evNum],['NitrogenFlow_' evNum]}
                        {'dillution', 'dillution'}
                        {'finalValveOnTime', 'fvOnTime'}
                        {'finalValveDuration', 'fvdur'}
                        {'finalValveTrigExh', 'fvtrig_on_exh'}
                        {'odor', 'odor'}
                        {'odorConcentration','odorconc'}
                        {'vial', 'vial'}
                        {'vialConcentration','vialconc'}];

    end %function extract_final_valve
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function evtFields= extract_laser_pars(evName)
        evNameSplit=strsplit(evName,'_');
        evNum = (evNameSplit{end});
        %generate the fields to get from the voyeur_table
        evtFields = [  
            {['laserIntensity_' evNum],['LaserIntensity_' evNum]}
            {['laserAmplitude_' evNum],['amplitude_' evNum]}
            {['pulseOnDur_' evNum],['pulseOnDur_' evNum]}
            {['pulseOffDur_' evNum],['pulseOffDur_' evNum]}
            {['pulseOnsetDelay_' evNum],['pulseOnDelay_' evNum]}
            {['trainLength_' evNum],['trainLength_' evNum]}
            {'laserMultiSniff', 'laser_multi_sniff'}
            {'laserTrigExh', 'laser_on_exh'}
            {'laserOnTime', 'laserontime'}
            {'numLasers','lasersnum'}
            ];
    end %function extract_final_valve
    
end

function table = h52m_table(table)
%gets an h5read table and re-format it to matlab-friendly data types
%very basic: just double for num classes
fieldsList = fields(table);
for jf = 1:numel(fieldsList)
    fieldName = fieldsList{jf};
    if isnumeric(table.(fieldName))
        % all numeric types go to double
        table.(fieldName) = double(table.(fieldName));
    elseif ischar(table.(fieldName))
        %the table goes to a cell array
        newCol = [];
        for il = 1:length(table.(fieldName))
            newCol = [newCol ; {deblank(table.(fieldName)(:,il)')}];
        end
        rmfield(table,fieldName);
        table.(fieldName) = newCol;
    end
end
end