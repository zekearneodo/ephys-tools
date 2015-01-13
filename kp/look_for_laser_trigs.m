function trigx = look_for_laser_trigs(mouse,sess,rec)
%finds laser shots of KPawakeM72_004, rec 'a'.
%cuts a chunk before and after, and creates a trial structure compatible
%with the one made by data_management_tools_032.

ft=file_tools();
bTrial=2000;
aTrial=2000;
tLaser=500; %time of the laser pulse after onset of trial;

if nargin > 1
    if isnumeric(mouse)
        vp.mouse = sprintf('%04d', mouse);
    end
    if isnumeric(sess)
        vp.sess = sprintf('%03d', sess);
    end
end

if nargin >2 && ~isempty(rec)
    if isnumeric(rec)
        vp.rec = sprintf('%02d', rec);
    end
end

fn=file_names(mouse,sess,rec);
q=load(fn.ss_sess_info);
irec=1;
thisRec=q.info.rec(irec);

sF=thisRec.sampling_freq;

fprintf('Getting the triggers of set of manual trials');
%open the file:
recDatFn=fullfile(fn.fold_ss_rec,sprintf('rec_%s.dat',rec))
dFid=fopen(recDatFn,'r');
trigChN=find(strcmpi({thisRec.chan.name},'Laser'));
%read the laser channel;
fprintf('++ Reading Laser channel from ephys file...\n');
fprintf('   ');
trigx=ft.read_analog_channel(dFid,trigChN,thisRec.nChan);
fclose(dFid);



if(mean(trigx)<0)
    thresh=0;
else
    thresh=mean(trigx)+0.5*range(trigx);
end
ntrigx=trigx;
ntrigx(trigx>thresh)=1;
ntrigx(trigx<thresh)=0;
%trigx=trigx-mean(trigx);
% figure;
% plot(trigx);
% hold on
%plot(ntrigx*range(trigx),'k')
trigOn=diff(ntrigx)>0;
trigOff=diff(ntrigx)<0;
% plot(diff(ntrigx)*range(trigx),'k')
% plot(trigOn*range(trigx),'r*')
% plot(trigOff*range(trigx),'go')
stampsOn=find(trigOn);
stampsOff=find(trigOff);
if(stampsOn(1)<ceil(sF*1)+1)
    stampsOn(1)=[];
    stampsOff(1)=[];
end
% if(stampsOff(end)>run.nSampl-ceil(Fs*1))
%     stampsOn(end)=[];
%     stampsOff(end)=[];
% end


% added by kp for trying to pull out laser amplitudes
amps = zeros(length(stampsOn),1);
for ia=1:length(stampsOn)
    
    amps(ia) = mean(trigx(stampsOn(ia)+1:stampsOff(ia)));
    
end
ampbins = min(amps)-50:10:max(amps)+50;
groups = histc(amps,ampbins);
figure; bar(groups)


stampsOn_time = round(stampsOn./sF.*1000);
load(fn.trial);
laserTimes = find([trial.laserTimes]);

ilas = 0; t_off = 0; 
las_used = []; stamps_notUsed = [];
for itri=1:length(trial)
    curr_trial = trial(itri);
    if curr_trial.odorConc==0 && curr_trial.laserTimes(1)>0
        
        closest_trig = find(abs(stampsOn_time - trial(itri).start)<3000);
        
        if isempty(closest_trig)
            sprintf('NO STAMP FOUND FOR TRIAL %i',itri)
            continue
        elseif closest_trig <= ilas
            error('WANTS TO USE A STAMP MORE THAN ONCE')
%         elseif closest_trig > (ilas+1+t_off)
%             stamps_notUsed = [stamps_notUsed (ilas+t_off+1):(closest_trig-1)]
%             t_off = t_off + closest_trig-ilas;
%             sprintf('no matching trials for trigs %i',stamps_notUsed)
        end
        
        las_used = [las_used closest_trig];
        laserON = stampsOn_time(closest_trig);
            NEW_laserOn = laserON-curr_trial.start;
        laserOFF = round(stampsOff(closest_trig)/sF*1000);
        
        trial_error = NEW_laserOn - curr_trial.laserTimes(1);
        if abs(trial_error) < 2500
            ilas = ilas+1;
            NEW_laserTimes(:,itri) = [NEW_laserOn; laserOFF-curr_trial.start];
            
            trial(itri).laserTimes = NEW_laserTimes(:,itri);
        else
            error('SHIT')
        end
        
        %% to check with neuroscope
        disp('---------------------------')
        itri
        trial_start = curr_trial.start
        NEW_laserOn

    else
        continue
    end
end






%that was the easy part. Now you've got the laser ticks.
%
trial=struct();
laserAmp=1050;
laserPower=1;
stimID=0;
run=1; %don't even care about the run right now.
odorTimes=[];
odorConc=0;
odorName='none';


for itrig=1:length(stampsOn)
    trial(itrig).rec = rec;
    trial(itrig).run = thisRec.run(find([thisRec.run.offset]<stampsOn(itrig),1)).num;
    trial(itrig).start      = round(stampsOn(itrig)/sF*1000)-tLaser;
    trial(itrig).duration   = aTrial+bTrial;
    trial(itrig).odorName   = odorName;
    trial(itrig).odorConc   = odorConc;
    trial(itrig).odorTimes  = odorTimes;
    trial(itrig).laserTimes = round([stampsOn(itrig);stampsOff(itrig)]/sF*1000)-trial(itrig).start;
    trial(itrig).laserAmp   = laserAmp;
    
    trial(itrig).laserPower = laserPower;
    trial(itrig).stimID     = stimID;
end
%save the trial structure
save(fn.trial, 'trial')


    %get_the_sniffs
    %gets the sniff crossings and dumps them into the trilas (as
    %sniffZeroTimes).
    %get the sniff
    dFid=fopen(recDatFn,'r');
    sniffChN=find(strcmpi({thisRec.chan.name},'Sniff'));
    %read one by one the ephys channels and compute the avg channel;
    fprintf('++ Reading Sniff channel from ephys file...\n');
    fprintf('   ');
    sniffx=ft.read_analog_channel(dFid,sniffChN,thisRec.nChan);
    fclose(dFid);
    
    %resample at 1Khz and find the positive-slope,zero crossings
    sniffxRes  = resample(sniffx,1000,sF);
    
    [bf, af] = butter(2, 20/500, 'low');
    
    sinffxRes  = filtfilt(bf, af, sniffxRes);
    sniffxRes  = -1*sniffxRes;
    sniffCross = find((diff(sniffxRes>0))>0);
    slopeCross = zeros(size(sniffCross));
    for iic=1:numel(sniffCross)
        ic=sniffCross(iic);
        slopeCross(iic)=max(sniffxRes(ic:min(ic+75,length(sniffxRes))))-min(sniffxRes(max(1,ic-75):ic));
    end
    slThresh=mean(findpeaks(sniffxRes,'minpeakheight',1e4))-.75*std(findpeaks(sniffxRes,'minpeakheight',1e4));
    keepCross=find(slopeCross>slThresh);
    sniffCross=sniffCross(keepCross);
    slopeCross=slopeCross(keepCross);
    downCross  = find((diff(sniffxRes<-40/2500*32768))>0);
    repeated = find(diff(sniffCross)<100);
    slopeCross(repeated)=[];
    sniffCross(repeated)=[];
    
    figure
    hold
    plot(sniffxRes,'r');
    plot(sniffCross,zeros(size(sniffCross)),'ob')
    %plot(downCross,-40/2500*32768*ones(size(downCross)),'ok')
    %plot(sniffCross(repeated),zeros(size(repeated)),'xg')
    stem(sniffCross,slopeCross,'xk')
    %plot(sniffCross(find(slopeCross>0)),zeros(size(find(slopeCross>0))),'k*');
    
    %now dump the sniffZerotimes in the trials
    %referred to the beginning of the trial
     for it=1:length(stampsOn)
        t0 = trial(it).start;
        t1 = t0 - bTrial;
        t2 = t0 + aTrial;
        
        
        in_sn = (sniffCross > t1)&(sniffCross < t2);
        

        
        trial(it).sniffZeroTimes      = sniffCross(in_sn) - t0;

        %find the pulseOffset (to check whether it makes some sense)
        tr=trial(it);
        trigSnif=max(tr.sniffZeroTimes(find(tr.sniffZeroTimes<tr.laserTimes(1,1))));
        if isempty(in_sn) || isempty(trial(it).sniffZeroTimes) || isempty(trigSnif)
            trial(it).trigSnif=-1
            trial(it).pulseOffset = -1;
            continue
        end
        trial(it).pulseOffset = tr.laserTimes(1,1) - trigSnif;
        trial(it).trigSnif = trigSnif;
    end
    save(fn.trial, 'trial')
   
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function Y = read_channel(fnData, ich, nch)
fidR = fopen(fnData, 'r');
fseek(fidR, 2*(ich-1), 'bof');
Y = fread(fidR, inf, 'int16', (nch-1)*2);
fclose(fidR);
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~