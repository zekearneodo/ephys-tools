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
inPar.addParamValue('samplingFrequency',19385,@(x) isscalar(x))

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