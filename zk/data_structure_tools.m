%Set of tools for extracting events and making trial structures.
%
%
function ds = data_structure_tools()
global ds

    ds.match_trial_numbers = @match_trial_numbers; %get trial numbers and match them with trial pin events
    ds.get_analog_events = @get_analog_events; %get occurrences of an analog event 
    ds.get_trial_numbers = @get_trial_numbers; % get the trial numbers
    ds.get_data_stream   = @get_data_stream; %get stream of data from binary file
    ds.get_info          = @get_info; %gets info structures for rec, run
    ds.check_list        = @check_list; %checks a list for repeated or skipped positions
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
%look for errors in the list
er = check_list(numStartPins);
if ~isempty(find(er==1, 1))
    warning('Some trial pins did not have a matching trial number');
end
if ~isempty(find(er==2, 1))
    warning('Some trial pins have several ');
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
    end
end

tn.table.trial_num   = [tn.events.value];
tn.table.trial_on    = [tn.events.on];
tn.table.trial_pin   = numStartPins;
tn.table.trial_start = tStartPins;
end

function get_fv(mouse,sess,rec,irun)
%gets the final valve time stamps
%gets the trial number that corresponds to that opening
%gets the odor parameters that correspond to that fv opening

    
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

function [trials, table] = get_voyeur_table(mouse,sess,rec,run)
%reads the table in a voyeur-generated behavior file for a run.
% Assumes there is only one table and that table is the trial data
% Doesn't care about the streams
% returns:
% trials : array of trial structures with data in mlab format (double and char)
% table  : h5-formatted table, as it was read.



[~, runInfo] = get_info(mouse,sess,rec,run);
fn = file_names(mouse,sess,rec,run);

fnam    = fullfile(fn.fold_rd_sess, runInfo.behav_data);
hinfo   = h5info(fnam);

n_tr = numel(hinfo.Groups);
table_name = hinfo.Datasets.Name;
table = h5read(fnam,['/', table_name]);
ff = fieldnames(table);
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