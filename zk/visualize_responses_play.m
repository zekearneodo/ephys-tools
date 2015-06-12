% Script for visualizing spikes aligned with the onset of stimuli
% based on data_managment_tools_03
% inits the parameters and the functions for analysis,
% leaves global vp with the parameters for all the functions to call
%has functions:
% -find_stimuli_set
% to play and get to a cell passport
function vr = visualize_responses_play(mouse,sess,rec,unit,sType)
global vp vd dm;  

%add the package of common files and folders to the path.
%here are the functions for file_naming, for instance, and any other
%functions that several script bundles will use.
%that is /baseFolder/ephysDataManagement/current/include
[~,computerName]=system('hostname');
if ~strcmp(strtrim(computerName),'flipper')
    includePath=fullfile(fileparts(pwd),'current','include');
    addpath(includePath);
end
% Functions of include that it uses:
%   - file_names(mouse,sess,rec,stat)

vr.visual_par_init     = @visual_par_init;
vr.make_stimuli_set    = @make_stimuli_set;
vr.view_rasters        = @view_rasters;
vr.view_lfp            = @view_lfp;
% pp.ss_prep        = @ss_prep;
% pp.resampling     = @resampling;

if nargin < 5
    sType='laser'
end

if nargin >1
dm = data_management_tools_031();
vp = visual_par_init(mouse,sess,rec,sType);
vd = make_stimuli_set(vp);

%if it goes withs spikes, do the raster visualization
%otherwies, we'll see.
if ~vp.surface && ~isempty(unit)
    view_rasters(sType,unit);
end
end

end

function vp=visual_par_init(mouse,sess,rec,sType)
% visual_par_init(mouse,sess,rec,[unit])
%initializes or resets the parameters for the visualization
%that will be stored in the global vp;
%makes the stimuli set
global dm vd;

vp.mouse = mouse;
vp.sess  = sess;
vp.rec   = rec;

%un    = 13;

%     fn.fold_sd_data  = fullfile(stat_disk, '');
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


vp.fn = file_names(mouse,sess,rec);
q=load(vp.fn.sess_info);
nrec=find(strcmp(rec,{q.info.rec.name}));

if strfind(q.info.rec(nrec).electrode_type,'surface');
    vp.surface=1;
else
    vp.surface=0;
end

if strcmp(sType,'laser')
    vp.t1=-200;
    vp.t2=400;
    vp.responseWindowDefault=100;
else
    vp.t1  = -200;
    vp.t2  =  400;
    vp.responseWindowDefault=200;
end
    
vp.nt  = vp.t2-vp.t1;
vp.bin = 20;
vp.t   = mean(reshape((vp.t1+1):vp.t2, vp.bin , vp.nt/vp.bin),1);
vp.smooth.wsize = round(vp.nt/8); % window size for gaussian smoothing of histo for plotting
vp.smooth.cutof = 100; % cutoff for the gausiann smoothing
vp.smooth.stdev = q.info.rec(nrec).sampling_freq/(2*pi*vp.smooth.cutof); % window size for gaussian smoothing of histo for plotting

q = load(vp.fn.trial);

%check if there is a trial correction function to be applied
trCFname=fullfile(vp.fn.fold_pr_sess,sprintf('%s_%s_trial_correct.m', vp.mouse,vp.sess));
if exist(trCFname)
    fprintf('Found trial correction function %s\n',trCFname);
    p = path;
    path(p,vp.fn.fold_pr_sess);
    fun=sprintf('%s_%s_trial_correct', vp.mouse,vp.sess);
    trCorrectFcn=eval(['@' fun ';'])
    path(p);
    trial=feval(trCorrectFcn,q.trial);
    vp.tr=trial;
else
    vp.tr = q.trial;
end

vp.stimTypes={'laser','odor'};
%parameters for determining uniqueness of stimulus
if isfield(vp.tr,'pulseOnsetDelay')
    vp.tr_type = 3;
    vp.par = {'odorName', 'odorConc', 'laserDur', 'laserPower', 'pulseOnsetDelay'}
elseif isfield(vp.tr,'pulseOffset')        % just KPawakeM72_004
    vp.tr_type = 2;
    vp.par = {'odorName', 'odorConc', 'laserDur', 'laserPower','pulseGroup','pulseOffset'}
else                                       % old sessions
    vp.tr_type = 1;
    vp.par = {'odorName', 'odorConc', 'laserDur', 'laserAmp'}
end

vd=make_stimuli_set(vp);
end

function vd=make_stimuli_set(vp)
% ====================================================================================
% finding stimulus set, and filtering by properties of stimulation
% eg: odorstim,lightstim

stim = struct('odorName', '', 'odorConc', 0, 'laserDur', 0, 'laserAmp', 0, 'in_tr', [])

for it=1:numel(vp.stimTypes)
    sType=vp.stimTypes{it}
    vd.(sType) = struct('list',[],'sort','')
end

% ====================================================================================
% creates set of stimuli to plot, separating by properties, depending on 
% eg: odorstim,lightstim
%TODO: turn this into an optional function that can be passed.


for it=1:numel(vp.stimTypes)
    sType=vp.stimTypes{it}
    vd.(sType) = struct('list',[],'sort','')
end

switch vp.tr_type
    case 1
        stim = struct('odorName', '', 'odorConc', 0, 'laserDur', 0, 'laserAmp', 0, 'in_tr', [],'odorInfo', {})
        tr=vp.tr;
    case 2
        stim = struct('odorName', '', 'odorConc', 0, 'laserDur', 0, 'laserPower', 0, 'in_tr', [],'odorInfo', {}, 'pulseOffset', [],'pulseGroup',[])
        tr=vp.tr;
    case 3
        stim = struct('odorName', '', 'odorConc', 0, 'laserDur', 0, 'laserPower', 0, 'in_tr', [],'odorInfo', {}, 'pulseOnsetDelay', [])
        tr=vp.tr;
end


if vp.tr_type==2                            % for KPawakeM72 sess_004 only
%         Group the laser pulse times according to sniff cycle
        ms_from_inh = [ 21 28;   68 75;   77 78;   101 108];  %these are the windows in ms!
        for jt = 1:numel(vp.tr)
            vp.tr(jt).pulseGroup = NaN;
        end
        used = [];
        for grp = 1:size(ms_from_inh,1)
            grp_ind = find([vp.tr.pulseOffset]>ms_from_inh(grp,1) & [vp.tr.pulseOffset]<ms_from_inh(grp,2));
            used = [used grp_ind];
            for jt = 1:numel(vp.tr)
                if ismember(jt,grp_ind)
                    vp.tr(jt).pulseGroup = grp;
                    vp.tr(jt).pulseOffset=round(mean(ms_from_inh(grp,:)));
                end
            end
        end
    tr=vp.tr(find([vp.tr.pulseGroup]>0)); 
end

%% Find each unique stimulus, and create list of trials corresponding to each

ns   = 1;
par=vp.par;
stim(1).odorName='';

for it = 1:numel(tr)
    if tr(it).start == 0
        continue
    end
        
    tr(it).laserDur = diff(tr(it).laserTimes);
    new = 1; ks = 0;
    while new&&(ks<ns)
        ks = ks + 1;
        same_stim = strcmp(tr(it).odorName, stim(ks).odorName);
        for ip = 2:numel(par)
            same_stim = same_stim&(tr(it).(par{ip}) == stim(ks).(par{ip}));
        end
        
        if same_stim
             new = 0;
             stim(ks).in_tr(end+1) = it;
        end
    end
    
    if new
        ns = ns + 1;
        for ip = 1:numel(par)
            stim(ns).(par{ip}) = tr(it).(par{ip});
        end
        %if the odorInfo structure wasn't made when trial structure was
        %made, make it now
        if isfield(stim,'odorInfo') && ~isfield(tr,'odorInfo')
            stim(ns).('odorInfo') = get_odor_info(tr(it).VoyeurParameters);
        end
        stim(ns).in_tr = it;
    end
end

stim = stim(2:end)
ns   = ns - 1;
%remove from the list the stimuli for which there are less than 10 trials
keepStim=[];
for is=1:numel(stim)
    if numel(stim(is).in_tr)>10
        keepStim=[keepStim is];
    end
end
stim=stim(keepStim);


%select stimuli types
for it=1:numel(vp.stimTypes)
    sType=vp.stimTypes{it};
    switch sType
        case 'odor'
            sSelect = find([stim.odorConc]>0);
            sSort   = 'odorName';
        case 'laser'
            sSelect = find([stim.odorConc]==0);
            if vp.tr_type==1
                sSort   = 'laserAmp';
            else
                sSort   = 'laserPower';
                %                 sSort   = 'pulseOnsetDelay';
            end
    end
    vd.(sType).list = sSelect;
    vd.(sType).sort = sSort;
end
vd.stim=stim;
end %function vd=make_stimuli_set()

function odorInfo=get_odor_info(voPar)
odorInfo=struct('name','','flows',[],'dillution',1,'vialConc',{},'vialId',{},'vial',{})

odorInfo(1).name      = voPar.odor;
odorInfo(1).flows     = [voPar.AirFlow_1 ;  voPar.NitrogenFlow_1];
odorInfo(1).dillution = [voPar.dillution];
odorInfo(1).vialConc  = [voPar.vialconc];
odorInfo(1).vialId    = [voPar.vial];

end

function view_rasters(sType,un)
global vd vp;

stim=vd.stim(vd.(sType).list);
sSort=vd.(sType).sort;

if isnumeric(stim(1).(sSort))
    [~, in_sort] = sort(-[stim.(sSort)]);
else
    [~, in_sort] = sort({stim.(sSort)});
end

stim = stim(in_sort);
orderedStim = [];
%if they are odors, sort by concentration
if strcmpi(sSort,'odorName')
    odorsList = unique({stim.odorName});
    odorsList(strcmpi({'not_an_event'},odorsList))=[];
    for io = 1:numel(odorsList)
        stimThisOdor = stim(strcmpi({stim.odorName},odorsList{io}))
        [~, in_sort] = sort(-[stimThisOdor.odorConc]);
        orderedStim = [orderedStim stimThisOdor(in_sort)];
    end
stim = orderedStim;    
end

ns=numel(stim);

for i=1:ns
    stim(i).lineStyle = '-';
    stim(i).lineColor = 'k';
    if stim(i).odorConc < 0.014
        stim(i).lineStyle = '-';
        stim(i).lineColor = 'r';
    end
end


rmax = 100;

%% =======================================================================
%   Plot raster and psth for defined units
%  =======================================================================

ifigRas=12;
ifigPsth=22;
gf = figure(ifigRas); clf
ip = 0;
tr=vp.tr;
t1=vp.t1;
t2=vp.t2;
bin=vp.bin;
t=vp.t;
smooth = vp.smooth;
nt=vp.nt;
tms = t1+1:t2;
nt=vp.nt;
fn=vp.fn;

% Get baseline rate for odor and laser trials             
%    averaging sime size window around inhalation onset 
%    for the 3 blank sniffs before each trial start
rateBase=zeros(1,nt);
nb=0; inhLenBase = []; inhLenStim = [];
sniffLengths   = zeros(numel(tr),3);
inhLengthsStim = zeros(numel(tr),4);
for kt = 1:numel(tr)
    trial=tr(kt);
    
%     if (kt==79)
%         kt
%     end
    if isempty(trial.sniffZeroTimes) || isempty(trial.spikeTimes)
        continue
    end
    sp = sort(round(vertcat(trial.spikeTimes{un})));

    nWindow=0;
    oneMoreWindow=1;
    while oneMoreWindow   % now averaging the [t1:t2] long windows that fit before the start of the trial.
        if nWindow==0
            %the first window is the one that ends when the firs sniff starts
            iFirstSniff = find(trial.sniffZeroTimes(1,:)>trial.odorTimes(1),1);
            t0Base = trial.sniffZeroTimes(1,iFirstSniff);

%             if trial.odorConc>0
%                 inhLenStim = [inhLenStim diff(trial.sniffZeroTimes(:,iFirstSniff))];
%                 %get the durations of the three whole sniffs after onset of
%                 %inhale
%                 sniffLengths(kt,:)   = diff(trial.sniffZeroTimes(1,iFirstSniff:iFirstSniff+3));
%                 inhLengthsStim(kt,:) = diff(trial.sniffZeroTimes(:,iFirstSniff:iFirstSniff+3));
%             end
%              % also get all the inhale durations
%             inhLenBase = [inhLenBase diff(trial.sniffZeroTimes(:,1:iFirstSniff-1))];
        else
            %for the subsequent windows it goes backwards
            if isempty(t0Base)
                break
            end
            iFirstSniff = find(trial.sniffZeroTimes(1,:)<t0Base-(vp.t2-vp.t1),1,'last');
            
            if isempty(iFirstSniff) || iFirstSniff<1
                break
            end
            
            t0Base = trial.sniffZeroTimes(1,iFirstSniff);
            spBase = sp((sp>t0Base+t1)&(sp<=t0Base+t2))-t0Base;
            rateBase(spBase-t1) = rateBase(spBase-t1)+1;
            nb=nb+1; %keep track of the number of segments added up to the baseline
            
            if t0Base-(vp.t2-vp.t1) < min(trial.sniffZeroTimes(1,:))
                oneMoreWindow=0;
            end
        end
        nWindow=nWindow+1;
    end
end
rateBase=rateBase/(nb); %divide by number of "trials"

%figure; hold on
%subplot(2,1,1)
%cdfplot(inhLenBase)
%xrange=0:300;
%title('no odor inhale duration')
%subplot(2,1,2)
%cdfplot(inhLenStim)
%xrange=0:300;
%title('first inhale duration')
%hold off

%if it is an odor trial, set the vp.responseWindow to the average duration
%of the first sniff


% spkBaseline = (1000/(0-t1))*sum(rateBase([t1:t2]<0))/nb;
%kp's accounting for spikes in the baseline

rateBaseHist = sum(reshape(rateBase, bin, nt/bin),1)/(bin/1000); %units is Hz
avgRateBase = mean(rateBase)*1000; %units is Hz



%% Now, for each stimulus find spikes in analysis window for each trial
%  Store data for raster (x,y) and histogram (rate, reshaped)
fg = figure();

for ks = 1:numel(stim)
    %if its odor trial the response window is the average duration of first
    %sniff
    if ~any(strcmpi(stim(ks).odorName,{'none','dummy'}))
        vp.responseWindow=mean(inhLenStim);
    else
        vp.responseWindow=responseWindowDefault;
    end
    timeObs=t1+1:t2;
    responseWindow=vp.responseWindow;
    spkCountBase = sum(rateBase(timeObs>0 & timeObs<responseWindow));
    nSpikesBase=sum(rateBaseHist(t>0 & t<responseWindow)*bin)/1000;
    nSpikesBaseAvg=avgRateBase*responseWindow/1000;
%     rateBase = smoothts(rateBase,'g',smooth.wsize,smooth.stdev);

    %     rmax=numel(stim(ks).in_tr);
    rate = zeros(1,vp.nt);
    
    rasters=[];
    x    = zeros(1e5,1);    y    = zeros(1e5,1);
    nsp  = 0;
    kt   = 0;
    
    if vp.tr_type == 1
        stim_intensity = stim(ks).laserAmp;
        stim_label = 'V';
    else
        stim_intensity = stim(ks).laserPower;
        stim_label = 'mW';
    end
    
%     upperTrials=min(20,numel(stim(ks).in_tr));
%     stim(ks).in_tr=stim(ks).in_tr(end-upperTrials+1:end);
    spkCount=[];
    for it = [stim(ks).in_tr]
        %get the spikeTimes for all the units
        sp = sort(round(vertcat(tr(it).spikeTimes{un})));
        
        %if its odor align with sniff, otherwise align with laser
        %presentation
        if stim(ks).odorConc >0 && strcmp(sType,'odor')
            %find first sniff after final valve open
            ii = find(tr(it).sniffZeroTimes(1,:)>tr(it).odorTimes(1)*1.009,1);
            if isempty(ii) || ii<1
                continue
            end
            t0 = tr(it).sniffZeroTimes(1,ii);
            stim_str = sprintf('%s \n %1.1e', tr(it).odorName(1:min(length(tr(it).odorName),10)), tr(it).odorConc);
            
        else
            t0 = tr(it).laserTimes(1);
            stim_str = sprintf('Laser: \n%2.1fV -%3.1fms', stim(ks).laserAmp/1000, stim(ks).laserDur);
        end
        kt = kt + 1;
        
        sp = sp((sp>t0+t1)&(sp<t0+t2))-t0;
        rate(sp-t1) = rate(sp-t1) + 1;
        n  = numel(sp);
        
        %another way of counting the number of spikes every trial
        %(how many spikes are from t0 to t0+200); but t0 is now the 0 in sp
        spkCount=[spkCount numel(sp((sp>0)&(sp<responseWindow)))];
        x(nsp+(1:n)) = sp(:);
        y(nsp+(1:n)) = kt*ones(1,n);
        nsp = nsp + n;
        
        % Calculate ISI for each sp in trial
        ISI(kt).ISI_raw = nan(1,nt);
        if n > 1
            ISI(kt).ISI_raw(1,sp(2:end)-t1) = diff(sp);  % ISI
            ISI(kt).ISI_raw(1,sp(1)-t1) = inf;
            ISI(kt).quickind = find(~isnan(ISI(kt).ISI_raw(1,:)));
            ISI(kt).quicklook = ISI(kt).ISI_raw(1,ISI(kt).quickind);
            
        else
            ISI(kt).quickind = [];
            ISI(kt).quicklook = [];
        end
        
    end
    x = x(1:nsp);
    y = y(1:nsp);
    rate = sum(reshape(rate, bin, nt/bin),1)/kt/(bin/1000);
    %     rateBase = sum(reshape(rateBase, bin, nt/bin),1)/kt*100;
    rate_pl = smoothts(rate,'g',smooth.wsize,smooth.stdev);
    
    
    baseline=mean(rate(t<0));
    pkThresh=baseline + 3.*std(rate(t<0));
    [pkRate,pkInd]=findpeaks(rate(t>0),'MINPEAKHEIGHT',pkThresh,'NPEAKS',1);
    pkTime=t(pkInd+find(t>0,1)-1);
    
    if isempty(pkTime)
        pkTime=0;
        pkRate=0;
    end
    
    nSpikes=(sum(rate(t>0 & t<responseWindow)*bin))/1000;
    
    % Caluculate response period from PSTH
    resp_thr = 0.7 * pkRate; % might have to be higher for less strong responses
      

% PLAY WITH DIFFERENT WAYS TO DEFINE THRESHOLD, SEE HOW MUCH LATENCY CHANGES

%         resp_pd = [min(t(rate > resp_thr)) max(t(rate > resp_thr))];
    isi_thr = 1/resp_thr*1000 ;% upper bound for ms since last spike
    isi_thr = 1/90*1000;
    
    % Pull stats on ISI from this response period
    ISItoplot = []; nsp = 0; 
    for jt = 1:kt
        underthr = find(ISI(jt).ISI_raw < isi_thr);
        ISI(jt).resp_spk_t = underthr;
        ISI(jt).resp_spk_t = underthr(underthr > 0-t1 & underthr < 400-t1);

        if ~isempty(ISI(jt).resp_spk_t)

            ISI(jt).respISI = ISI(jt).ISI_raw(1,underthr);
            ISItoplot = [ISItoplot ISI(jt).respISI];
            ISI(jt).nSpks = numel(ISI(jt).resp_spk_t); % note this does NOT count the spike that precedes the first under-threshold ISI
            ISI(jt).avgISI = mean(ISI(jt).ISI_raw(1,underthr));
            
            ISI(jt).first = ISI(jt).resp_spk_t(1) + t1;
            
        else
            ISI(jt).respISI = [];
            ISI(jt).nSpks = 0; % note this does NOT count the spike that precedes the first under-threshold ISI
            ISI(jt).avgISI = nan;
            ISI(jt).first = nan;
        end
        
        % Save raster info for spikes used in ISI measures
        nn = numel([ISI(jt).resp_spk_t]); 
        xx(nsp+(1:nn)) = [ISI(jt).resp_spk_t] +t1;   %%% histo of xx maybe related to pref sniff phase?? check with sniff stats, and with phase-triggered stimulus presentations
        yy(nsp+(1:nn)) = jt*ones(1,nn);
        nsp = nsp + nn;
    end
    
    ISI(:).first;
    pkTime = pkTime;
    latency = mean([ISI(isfinite([ISI(:).first])).first]);
    jitter = std([ISI(isfinite([ISI(:).first])).first]);
    avg_nSpks = mean([ISI(:).nSpks]);
    avg_ISI = mean([ISI(isfinite([ISI(:).avgISI])).avgISI]);
    
    vd.resp(ks).rateBase        = rateBase;
    vd.resp(ks).spkCountBase    = spkCountBase;
    vd.resp(ks).nSpikesBase     = nSpikesBase;
    vd.resp(ks).responseWindow  = responseWindow;
    %
    vd.resp(ks).ntrial          = numel(stim(ks).in_tr);
    vd.resp(ks).spkCount        = spkCount;
    vd.resp(ks).nSpikes         = nSpikes;
    vd.resp(ks).nSpk            = mean(spkCount);
    vd.resp(ks).nSpkErr         = std(spkCount);
    vd.resp(ks).latency         = pkTime;
    vd.resp(ks).maxFR           = pkRate;
    vd.resp(ks).latencyISI      = latency;
    vd.resp(ks).jitterISI       = jitter;
    vd.resp(ks).nSpksISI        = avg_nSpks;
    vd.resp(ks).respavgISI      = avg_ISI; 
    vd.resp(ks).vp              = vp;
    vd.resp(ks).inhLenghts      = inhLengthsStim;
    vd.resp(ks).inhLenghtsBase  = inhLenBase;
    
    vd.resp(ks).stim            = stim(ks);
    % variables and parameters for plotting raster
    vd.resp(ks).rateBaseHist    = rateBaseHist;
    vd.resp(ks).raster.x        = x;
    vd.resp(ks).raster.y        = y;
    vd.resp(ks).rate            = rate;
    vd.resp(ks).t               = t;
    vd.resp(ks).fg              = fg;
    vd.resp(ks).sType           = sType;
    %%%
    
    
    %plots this psth
    ip = ip + 1;
    gf = figure(ifigRas)
    gs = subplot(2,numel(stim), ip);
    yt  = max([1, floor(kt/10)*10]);
    plot(x,y, '.', 'MarkerSize',7), hold on
    plot([0,0], [0, kt+1], '--k'), hold off
    set(gs, 'XLim', [t1, t2], ...
        'YLim',      [0, kt+1], ...
        'YTick',     [0, yt], ...
        'XTick',     [-200,0,200], ...
        'FontSize',  10);
    title(stim_str, 'FontSize', 10, 'FontWeight', 'bold')
    
    
    gs = subplot(2,numel(stim), ip+numel(stim));
    plot(t, rate, ...
        'LineWidth', 1, ...
        'Color',     stim(ks).lineColor, ...
        'LineStyle', stim(ks).lineStyle),
    hold on
    %plot(t,ones(size(t))*baseline,'.k')
    %plot(t,ones(size(t))*pkThresh,'b')
    if strcmp(sType,'odor')
        plot(t,rateBaseHist,'LineStyle','--','Color',[.5,0.5,0.5])
    end
    plot(pkTime,pkRate,'k*')
    plot([0,0], [0, rmax], '--k'), hold off
    set(gs, 'XLim',     [t1, t2], ...
        'XTick',         [-200,0,200], ...
        'YLim',          [0, rmax], ...
        'YTick',         0:20:rmax, ...
        'FontSize',      10);
    
    %plots of psths altogether
    if ks>1
        figure(ifigPsth)
        gs0 = subplot(1,1,1);
        plot(t, rate, ...
            'LineWidth', 2, ...
            'Color',     stim(ks).lineColor, ...
            'LineStyle', stim(ks).lineStyle),
        hold on
    end
    
    
end

gf = figure(ifigRas)
sUnits=num2str(un(1));
for iu=2:numel(un)
    sUnits=sprintf('%s,%2d',sUnits,un(iu));
end
suptitle(sprintf('Mouse %s - session %3d - rec %s units %s',num2str(vp.mouse),str2num(vp.sess),vp.rec,sUnits));
%saves the response data (the plot and the response structure)
if ~exist(fn.fold_an_mouse, 'dir')
    mkdir(fn.fold_an_mouse)
end
if ~exist(fn.fold_an_sess, 'dir')
    mkdir(fn.fold_an_sess)
end

fileBase=[fn.basename_an sType '_units'];
for iu=1:numel(un)
    fileBase=[fileBase num2str(un(iu),'%02d')]
end

save(fullfile(fn.fold_an_sess,sprintf('%s_resp.mat',fileBase)), '-struct', 'vd','resp')
print(gf,'-depsc',fullfile(fullfile(fn.fold_an_sess,sprintf('%s_resp.eps',fileBase))))


if ks>1
    figure(ifigPsth)
    gs0 = subplot(1,1,1);
    subplot(gs0)
    plot([0,0], [0, rmax], '--k'), hold off
    set(gs0, 'XLim',     [t1, t2], ...
        'XTick',         -200:100:200, ...
        'YLim',          [0, rmax], ...
        'YTick',         0:20:rmax, ...
        'FontSize',      10);
end
end

function view_lfp(sType,chan,trialsPlot)
%views the trace of one channel, for firstTrials trials.
%if chan is a string, it is the channel name
%if chan is a number, it is the number of channel in the fn.rsm_data struct
%array.

global vd vp;
close all;
ip = 0;
tr=vp.tr;
t1=-300;
t2=500;
nt=vp.nt;
fn=vp.fn;

if ~exist('trialsPlot','var') || isempty('firstTrials')
    trialsPlot=1:9;
    numTrials=numel(trialsPlot);
end

tPlot=zeros(trialsPlot,2);
tVec=(t1:t2);

%get the list of channels
chanList=whos('-file',vp.fn.rsm_data);

%decide which channel to work on
if isnumeric(chan)
    %it is the number of channel
    chanNumber=chan;
elseif ~isempty(str2num(chan))
    %it is the number of channel but in string format
    chanNumber=str2num(chan);
else
    %it is the name of the channel
    chanNumber=find(strcmpi(chan,{chanList.name}));
end
   
%get the channel
theChan=load(vp.fn.rsm_data,'-mat',chanList(chanNumber).name);
chanData=theChan.(chanList(chanNumber).name);

%get the sniff
% snifNumber=find(strcmpi('sniff',{chanList.name}));
% theChan=load(vp.fn.rsm_data,'-mat',chanList(snifNumber).name);
% snifData=theChan.(chanList(snifNumber).name);

stim=vd.stim(vd.(sType).list);
ns=numel(stim);
sSort=vd.(sType).sort;

if isnumeric(stim(1).(sSort))
    [~, in_sort] = sort(-[stim.(sSort)]);
else
    [~, in_sort] = sort({stim.(sSort)});
end
stim = stim(in_sort)

ns=numel(stim);

for i=1:ns
    stim(i).lineColor = 'r';
end

for ks = 1:numel(stim)
    kt   = 0;
    
    ifigLfp=figure(ks);
    gf = figure(ks); clf
    stim_str = sprintf('%s %1.1e - Ch %s', ...
        tr(stim(ks).in_tr(1)).odorName(1:min(length(tr(stim(ks).in_tr(1)).odorName),10)),...
        tr(stim(ks).in_tr(1)).odorConc,...
        chanList(chanNumber).name);
    
    
    trialsToPlot=trialsPlot(trialsPlot<=numel(stim(ks).in_tr));
    numTrials=numel(stim(ks).in_tr(trialsToPlot));
    cols=ceil(sqrt(numTrials));
    isp=0;
    
    for it = stim(ks).in_tr(trialsToPlot)
        %if its odor align with sniff, otherwise align with laser
        %presentation
        if stim(ks).odorConc >0 & strcmp(sType,'odor')
            %find first sniff after final valve open
            ii = find(tr(it).sniffZeroTimes(1,:)>tr(it).odorTimes(1),1);
            if isempty(ii) | ii<1
                continue
            end
            t0 = tr(it).sniffZeroTimes(1,ii);
            
        else
            t0 = tr(it).laserTimes(1);
            stim_str = sprintf('Laser: \n%2.1fV -%3.1fms', stim(ks).laserAmp/1000, stim(ks).laserDur);
        end
        kt = kt + 1;
        
        
%         sp = sp(t0+t1)-t0;
        tPlot(it,:)=[t1 t2]+t0*0+tr(it).start;
        
        %plots this chunk
        isp=isp+1;
        gs = subplot(cols,cols,isp);
        yt  = max([1, floor(kt/10)*10]);
        plot(tVec,chanData(tPlot(it,1):tPlot(it,2)), 'MarkerSize',7), hold on
        plot([0,0], [0, kt+1], '--k'), hold off
        set(gs, 'XLim', [t1, t2], ...
            'XTick',     [-200,0,200], ...
            'FontSize',  10);
    end
        suptitle(stim_str)
end

end