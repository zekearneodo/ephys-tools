%Feb 2014.
% This is the script on which the trial_shift_estim of
% data_management_tols_031 is based.
%It is meant to debug when the trial alignment looks funny.
%It is also the lab to test newer features, like handling skipped trials
%and "negative" shift.
%As it is now, the shift (positive) means that beh. trials were not
%recorded; so the first iShift trials of the beh are dropped before trying
%the alignment, drift estimate and error calculation.
%ToDo:
% -"Negative" shift
% - Skip trials
% 
clear all
close all

dm=data_management_tools_032;
recs={'a','b','c','d','e','f','g'};

irec=7;
irun=3;
iShift=6;


eTrig_chan = 'TrPin';
bTrig_chan = 'runTrialStart';

fn=dm.file_names('ZKawakeM72',001,recs{irec});

%load info
q = load(fn.sess_info);
info = q.info;

%load trial structure (for the rec)
q = load(fn.trial);
trial = q.trial;

%load ephys sync channel (for the rec)
q = load(fn.rsm_data, eTrig_chan);

%load behavior trig channel (for the run; timing relative to the run)
if ~ strcmp('runTrialStart',bTrig_chan)
    tt = [1;1]*[trial.runTrialStart] + [trial.(bTrig_chan)];
    bSync = tt(1,:);
else
    bSync=[trial.runTrialStart];
end

% run-specific

%get the ephys channel for this run
run = info.rec(irec).run(irun);
run_start = run.start;
run_dur   = run.duration;

eTRun=q.(eTrig_chan)(run.start+(1:run.duration-1));
eTRun(1:1.2E5)=0;
%trial beginnings
eThresh=min(eTRun)+range(eTRun)*1/2
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

figure
plot(eTRun)
hold
plot(eSyncRunStart,ones(size(eSyncRunStart))*15000,'r*')
plot(eSyncRunEnd,ones(size(eSyncRunEnd))*15000,'g*')
eSyncRun=[eSyncRunStart;eSyncRunEnd];
plot(eSyncRun,ones(size(eSyncRun))*15000);

%get the good trials for this run
in=find( ([trial.run] == irun) & ([trial.runTrialDur]>0));
bSyncRun = [bSync(in);bSync(in)+[trial(in).runTrialDur]];
%discard the trials with null duration
% inGood   = [trial(in).runTrialDur]>0;
% bSyncRun=bSyncRun(:,inGood)

bSyncRunPlot=bSyncRun-(bSyncRun(1,1+iShift)-eSyncRun(1,1));
nn=min(length(eSyncRun),length(bSyncRun));
alfa=eSyncRun(1,1:nn)/bSyncRunPlot(1,1:nn);
plot(bSyncRunPlot,ones(size(bSyncRunPlot))*10000,'*')


%align them for one shift

if (iShift>=0) %positive shift is ephys started after behavior (first beh trials skipped)
bSR=bSyncRun(:,1+iShift:end);
in=in(1+iShift:end);
eSR=eSyncRun;
elseif iShift<0 %negative shift is from ephys end (
eSR=eSyncRun(:,end-(length(eSyncRun)-iShift):end)
end

%set to beginning of train
eT0=eSR(1,1);
bT0=bSR(1,1);

eSR=eSR-eT0;
bSR=bSR-bT0;

nn=min(length(eSR),length(bSR));
eSR=eSR(:,1:nn);
bSR=bSR(:,1:nn);
ind=1:nn;

trRun=in(ind)
alpha=eSR(1,:)/bSR(1,:);
eta=sqrt(norm(eSR-alpha*bSR))
errorVec=sum(eSR-alpha*bSR,1);

figure
plot(eTRun(eT0:end))
hold
plot(eSR,ones(size(eSR))*15000)
plot(bSR*alpha,ones(size(eSR))*10000,'+')



figure
plot(eSR,bSR*alpha,'*')


