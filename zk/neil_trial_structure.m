function trials=neil_trial_structure(mouse,sess,rec)
%gets a trial structure for a mouse,sess,rec and writes it in a format
%suitable for neil
fn=file_names(mouse,sess,rec);
load(fn.trial);
load(fn.sess_info);

recInfo =  info.rec(strcmpi(rec,{info.rec.name}));
rsmSniff = load(fn.rsm_data,'Sniff');

%pick only the odor trials
trial(strcmpi({'none'},{trial.odorName}))=[];
trial(strcmpi({'empty'},{trial.odorName}))=[];
t1 = -200;
t2 = 2500;

trials = [];
for it=1:numel(trial)
    ttr = trial(it);
    tr.start      = ttr.start+t1;
    if isempty(tr.start) || ttr.odorTimes(1)<=0
        continue
    end
    tVec          = ttr.start + (t1:t2-1); %absolute times of the trial chunk
    tr.id         = sprintf('%s_%03d_%s_%d',mouse,sess,rec,ttr.start);
    tr.odorName   = ttr.odorName;
    tr.odorConc   = ttr.odorConc;
    tr.flow       = rsmSniff.Sniff(tVec);
    %the stimulus vector
    tr.stim       = zeros(1,t2-t1);
    tr.odorTimes  = ttr.odorTimes - t1 ;
    stimTimes  = tr.odorTimes;
    stimTimes([tr.odorTimes]>t2-t1)=t2-t1;
    stimTimes([tr.odorTimes]<0)=0;
    tr.stim(stimTimes(1):stimTimes(2))=1;

    %the sniff zero times
%     if it==105
%     disp(it)
%     end
    %some debugging:
    %plot the sniff, the stimulus and the sniffZero times for this trial
    %and check how they align:
%     figure
%     hold on
%     plot(tr.flow*-1)
%     plot(tr.stim*10000,'r')
%     plot(ttr.sniffZeroTimes(1,:)-t1,zeros(size(ttr.sniffZeroTimes(1,:))),'k*')
%     plot(ttr.sniffZeroTimes(2,:)-t1,zeros(size(ttr.sniffZeroTimes(2,:))),'go')
%     plot(ttr.sniffParabZeroTimes(2,:)-t1,zeros(size(ttr.sniffParabZeroTimes(2,:))),'mx')
    sn = zeros(size(tr.stim));
    try
    if ~any(isnan(ttr.sniffParabZeroTimes(2,:)))
        spZeros=[];
        spZeros(2,:) = round((ttr.sniffParabZeroTimes(2,:)));
        spZeros(1,:) = round((ttr.sniffZeroTimes(1,:)));
        tr.sniffZeroTimes = ttr.sniffZeroTimes - t1;
        sz = spZeros - t1*ones(size(spZeros));
        firstSnif = find(sz>0,1,'first');
        lastSnif  = find(sz<t2,1,'last');
        sn(1:sz(firstSnif))=bitget(firstSnif-1,1);
        for is=firstSnif:lastSnif-1
            sn(sz(is):sz(is+1))=bitget(is,1);
        end
        sn(sz(lastSnif):end)=bitget(lastSnif,1);
        sn(sn==0)=-1;
        tr.sniffPhase = sn;
        trials = [trials tr];
    end
    catch
        warning('Error in sniff phases of trial %d - %d',ttr.run,ttr.runTrialNum);
        tr.sniffPhase = nan;
    end
end

end
