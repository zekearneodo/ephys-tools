%Set of tools for extracting events and making trial structures.
%
%
function ds = data_structure_tools()
global ds

    ds.get_analog_events = @get_analog_events;
    ds.get_data_stream   = @get_data_stream;


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

function get_trials(mouse,sess,rec,irun)
%gets the time stamps of trial pins (on/off)
%gets the trial number that corresponds to that pin
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
inPar.addParamValue('figures','plot',@(x) ischar(x) && any(strcmpi(x,{'plot','noplot'})));
inPar.addParamValue('aqSystem','JFRC',@(x) ischar(x) && any(strcmpi(x,{'JFRC','intan'})));
inPar.addParamValue('threshold',0.25,@(x) isscalar(x));
inPar.addParamValue('chanGain',1,@(x) isscalar(x));

inPar.parse(chanName,mouse,sess, rec,irun,varargin{:});
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

function stream = get_data_stream(chanName, mouse,sess,rec,irun,varargin)
%read a stream of data from a binary file
        %read a channel with recorded events
        %eventChan = structure with channel info
        %stream    = structure with fields 
        %   data = the raw data in the channel
        %   int2Volt = conversion factor from integer value to voltage
        %   threshold = threshold (for on/off extraction) in integers
        [recInfo,runInfo] = get_info(mouse,sess,rec,irun);
        
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