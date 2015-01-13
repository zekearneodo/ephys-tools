function [onEvents]=get_analog_events(chanName,mouse, sess, rec,irun,varargin)
%tunrs an analog channel with events (Hi/Lo) into a semi-digital:
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