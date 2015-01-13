function laserEvents=get_laser(mouse, sess, rec,irun,varargin)
%gets the laser events (higher than 0.75V) out of the raw data file of a
%run


inPar=inputParser;
inPar.addRequired('mouse')
inPar.addRequired('sess')
inPar.addRequired('rec',@ischar)
inPar.addRequired('irun',@isscalar)
inPar.addParamValue('figures','plot',@(x) ischar(x) && any(strcmpi(x,{'plot','noplot'})));
inPar.addParamValue('aqSystem','JFRC',@(x) ischar(x) && any(strcmpi(x,{'JFRC','intan'})));
inPar.addParamValue('chanName','Laser',@(x) ischar(x));
inPar.addParamValue('chanGain',1,@(x) isscalar(x));

inPar.parse(mouse,sess, rec,irun,varargin{:});
figsOn=inPar.Results.figures;
laserChanName=inPar.Results.chanName;

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
laserChan=recInfo.chan(strcmpi(laserChanName,{recInfo.chan.name}));
fprintf('Reading laser channel (%s)...',laserChan.name);
runFID=fopen(runRawFile,'r');
laserData=read_analog_channel(runFID,laserChan.num_rd,recInfo.nChan_rd);
fclose(runFID);
fprintf('done.\n');

%it is easier to work with integers, so get the Voltage values to integers
volt2int=32768/(recInfo.rangeMax); % in V
int2volt=1./volt2int;

%the minimum detectable peak is 0.5v
minPeak=ceil(0.1*volt2int);
laserOnOff=zeros(size(laserData));
laserOnOff(laserData>0.1*volt2int)=1;

laserOnSample  = find(diff(laserOnOff)==1);
laserOffSample = find(diff(laserOnOff)==-1);

if laserOffSample(1)<laserOnSample(1)
    laserOffSample(1)=[];
end
if laserOnSample(end)>laserOffSample(end)
    laserOnSample(end)=[];
end

laserSample=[laserOnSample';laserOffSample'];

%get the mean values of every laser on segment
laserAmp   =arrayfun(@(x,y) mean(laserData(x+1:y-1)),laserSample(1,:),laserSample(2,:));

if figsOn
    figure
    hold on
    plot(laserData)
    %plot(laserOnOff*40000,'m')
    plot(laserOnSample,laserAmp,'r*')
    plot(laserOffSample,laserAmp,'go')
end

laserEvents=[laserSample;laserAmp*int2volt/inPar.Results.chanGain];


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function Y = read_analog_channel(fid, ich, nch)
        fseek(fid, 2*(ich-1), 'bof');
        Y = fread(fid, inf, 'int16', (nch-1)*2);
    end

end