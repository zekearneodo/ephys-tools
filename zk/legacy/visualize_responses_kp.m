%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script for visualizing spikes aligned with the onset of stimuli

% outputs vr and response struct
% global vp contains the parameters for all the functions to call

%has functions:
% -visual_par_init
% -make_stimuli_set
% -view_rasters (default)
% -view_lfp (not called by default)

function [vr, resp] = visualize_responses_kp(mouse,sess,rec,unit,sType)
global vp vd;  

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

if nargin < 5
    sType='laser'
end

vp = visual_par_init(mouse,sess,rec,sType);
vd = make_stimuli_set(vp);

if ~vp.surface && ~isempty(unit)
    resp = view_rasters(sType,unit);  %if spike data collected, do the raster visualization
end
end

function vp=visual_par_init(mouse,sess,rec,sType)
% visual_par_init(mouse,sess,rec,[unit])
%initializes or resets the parameters for the visualization
%that will be stored in the global vp;
%makes the stimuli set
global vd;

vp.mouse = mouse;
vp.sess  = sess;
vp.rec   = rec;

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
vp.bin = 10;
vp.t   = mean(reshape((vp.t1+1):vp.t2, vp.bin , vp.nt/vp.bin),1);
vp.smooth.wsize = round(vp.nt/8);   %window size for gaussian smoothing of histo for plotting
vp.smooth.cutoff = 100;   %cutoff for the gaussian smoothing
vp.smooth.stdev = q.info.rec(nrec).sampling_freq/(2*pi*vp.smooth.cutoff); %std for gaussian smoothing of histo for plotting

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

vd=make_stimuli_set(vp);
end %vp=visual_par_init(mouse,sess,rec,sType)

function vd=make_stimuli_set(vp)
% ====================================================================================
% Creates the set of unique stimuli to plot.
% - Based on fields of trial struct, sets up fields for stim struct
% - Finds all unique stimuli and creates lists of corresponding trials
% - 

%TODO: turn this into an optional function that can be passed.


%%%%%%%%%%%%
%%% Based on the fields of the trial structure, define appropriate stimulus
%%% parameters and create stim structure.
tr = vp.tr;

vp.stimTypes={'laser','odor'};
if isfield(vp.tr,'pulseOnsetDelay')
    vp.tr_type = 3;
    vp.par = {'odorName', 'odorConc', 'laserDur', 'laserPower', 'pulseOnsetDelay'}
elseif isfield(vp.tr,'pulseOffset')        
    vp.tr_type = 2;                        % just KPawakeM72_004
    vp.par = {'odorName', 'odorConc', 'laserDur', 'laserPower', 'pulseGroup','pulseOffset'}
else               
    vp.tr_type = 1;                        % old sessions
    vp.par = {'odorName', 'odorConc', 'laserDur', 'laserAmp'}
end

for it=1:numel(vp.stimTypes)
    sType=vp.stimTypes{it};
    vd.(sType) = struct('list',[],'sort','');
end

switch vp.tr_type
    case 1
        stim = struct('odorName', '', 'odorConc', 0, 'laserDur', 0, 'laserAmp', 0, 'in_tr', [], 'stim_str', '','odorInfo', {});
    case 2
        stim = struct('odorName', '', 'odorConc', 0, 'laserDur', 0, 'laserPower', 0, 'in_tr', [], 'pulseOffset', [],'pulseGroup',[], 'stim_str', '','odorInfo', {});
    case 3
        stim = struct('odorName', '', 'odorConc', 0, 'laserDur', 0, 'laserPower', 0, 'in_tr', [], 'pulseOnsetDelay', [], 'stim_str', '','odorInfo', {});
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

%%%%%%%%%%%%
%%% Find each unique stimulus and make list of trials corresponding to each

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

%remove the stimuli for which there are less than 10 trials
keepStim=[];
for is=1:numel(stim)
    if numel(stim(is).in_tr)>10
        keepStim=[keepStim is];
    end
end
stim=stim(keepStim);

%%%%%%%%%%%%
%%% Select stimuli and order for analysis, depending on sType argument
for it=1:numel(vp.stimTypes)
    sType=vp.stimTypes{it};
    switch sType
        case 'odor'
            sSelect = find([stim.odorConc]>0);
            sSort   = 'odorName';
            for is=1:numel(stim)
                stim_str = sprintf('%s \n [%1.1e]', stim(is).odorName, stim(is).odorConc);
                stim(is).stim_str = stim_str;
            end
            
        case 'laser'
            sSelect = find([stim.odorConc]==0);
            if vp.tr_type==1
                sSort   = 'laserAmp';
                for is=1:numel(stim)
                    stim_str = sprintf('Laser: \n%2.1fV, %3.1fms', stim(is).laserAmp, stim(is).laserDur);
                    stim(is).stim_str = stim_str;
                end
            else
                sSort   = 'laserPower';
                %                 sSort   = 'pulseOnsetDelay';
                for is=1:numel(stim)
                    stim_str = sprintf('Laser: \n%2.1fmW, %3.1fms', stim(is).laserPower, stim(is).laserDur);
                    stim(is).stim_str = stim_str;
                end
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

function resp = view_rasters(sType,un)
global vd vp;

stim=vd.stim(vd.(sType).list);
sSort=vd.(sType).sort;

if isnumeric(stim(1).(sSort))
    [~, in_sort] = sort(-[stim.(sSort)]);
else
    [~, in_sort] = sort({stim.(sSort)}); % here is where one could sort by affinity
end

stim = stim(in_sort);
ns=numel(stim);

for i=1:ns
    stim(i).lineStyle = '-';
    stim(i).lineColor = 'k';
    if stim(i).odorConc < 0.014
        stim(i).lineStyle = '-';
        stim(i).lineColor = 'r';
    end
end




%% =======================================================================
%   Plot raster and psth for requested units
%  =======================================================================

ifigRas=12;
ifigPsth=22;
gf  = figure(ifigRas); clf
ip  = 0;
tr  =vp.tr;
t1  =vp.t1;
t2  =vp.t2;
tms = t1+1:t2;  %time from t1:t2 in ms, not reshaped
t   =vp.t;      %time from t:t2, reshaped according to bin size
bin =vp.bin;
smooth = vp.smooth;
nt  =vp.nt;
fn  =vp.fn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            Get baseline rate              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Averaging spikes in same size window (t1:t2) around inhalation onsets 
% for a few blank sniffs before each trial start.

rateBase=zeros(1,nt);
nb=0; inhLenBase = []; inhLenStim = []; timeToFirstInhOdor = []; laserTimeAfterInhDetected = [];
nSniffs=4;
sniffLengths   = nan(numel(tr),nSniffs);
inhLengthsStim = nan(numel(tr),nSniffs+1);
for kt = 1:numel(tr)
    trial=tr(kt);
    if isempty(trial.sniffZeroTimes) || isempty(trial.spikeTimes)
        continue 
    end
    sp = sort(round(vertcat(trial.spikeTimes{un})));

    nWindow=0;
    oneMoreWindow=1;
    while oneMoreWindow 
        % This window is the first that that is part of the stimulus period
        if nWindow==0   
            if trial.odorConc>0
                iFirstSniff = find(trial.sniffZeroTimes(1,:)>trial.odorTimes(1),1,'first'); 
                t0Base = trial.sniffZeroTimes(1,iFirstSniff);  %time of first sniff after stimulus
                timeToFirstInhOdor = [timeToFirstInhOdor t0Base-trial.odorTimes(1)];
                
            elseif trial.laserAmp>0 && strcmp(trial.odorName,'none')
                iFirstSniff = find(trial.sniffZeroTimes(1,:)<trial.laserTimes(1),1,'last'); %t0 is TRIAL start (not FV or laser), could correct based on stim type
                t0Base = trial.sniffZeroTimes(1,iFirstSniff);  %time of first sniff after stimulus
                laserTimeAfterInhDetected = [laserTimeAfterInhDetected trial.laserTimes(1)-t0Base];
            end
            % Respiration statistics:
                %get duration of first inh after stim onset
            inhLenStim = [inhLenStim diff(trial.sniffZeroTimes(:,iFirstSniff))];
                %get durations of the three sniff cycles and inhalations
                %after stimulus onset
                for nS=nSniffs-(0:nSniffs-1)
                    try
                        sniffLengths(kt,1:nS) = diff(trial.sniffZeroTimes(1,iFirstSniff:iFirstSniff+nS));
                        inhLengthsStim(kt,1:nS+1) = diff(trial.sniffZeroTimes(:,iFirstSniff:iFirstSniff+nS));
                        break
                    catch
                        continue
                    end
                end
                %get all pre-stimulus inhale durations
            inhLenBase = [inhLenBase diff(trial.sniffZeroTimes(:,1:iFirstSniff-1))];
        
        % Now add up the spikes in as many non-overlapping windows as fit in the pre-stimulus period
        else
            iFirstSniff = find(trial.sniffZeroTimes(1,:)<t0Base-nt,1,'last');
            
            if isempty(iFirstSniff) || iFirstSniff<1
                break
            elseif trial.sniffZeroTimes(2,iFirstSniff)>0
                break
            end
            
            t0Base = trial.sniffZeroTimes(1,iFirstSniff);
            spBase = sp((sp>t0Base+t1)&(sp<=t0Base+t2))-t0Base;
            rateBase(spBase-t1) = rateBase(spBase-t1)+1;
            nb=nb+1; %keep track of the number of segments added to baseline
            
            if t0Base-nt < min(trial.sniffZeroTimes(1,:))
                oneMoreWindow=0;
            end
        end
        nWindow=nWindow+1;
    end
end
rateBase=rateBase/nb; %divide by number of segments used 

%% DEBUG THIS SOON!!
laseroffset=laserTimeAfterInhDetected 
%%

% plot respiration statistics
figure; hold on
subplot(3,1,1)
cdfplot(inhLenBase)
xrange=0:300;
title('no odor inhale duration')
subplot(3,1,2)
cdfplot(inhLenStim)
xrange=0:300;
title('first inhale with stim duration')
subplot(3,1,3)
cdfplot(timeToFirstInhOdor)
xrange=0:600;
title('latency to first inh after FV on')
hold off


    

%% =======================================================================
%%  For each stimulus find spikes in analysis window for each trial
%   Store data for raster (x,y) and histogram (rate, reshaped)

for ks = 1:numel(stim)
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx = sprintf('Stimulus #%i is %s',ks, stim(ks).stim_str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    First, calculate baseline properties   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rateBaseHist = sum(reshape(rateBase, bin, nt/bin),1)/(bin/1000); %units is Hz
    % rateBaseHist_smooth = smoothts(rateBaseHist,'g',smooth.wsize,smooth.stdev);
avgRateBase = mean(rateBase)*1000 %units is Hz

%if stim is an odor, set the response window to the average duration of first sniff
if ~any(strcmpi(stim(ks).odorName,{'none','dummy'}))
    vp.responseWindow=mean(inhLenStim);
else
    vp.responseWindow=responseWindowDefault;
end
responseWindow = vp.responseWindow;

spkCountBase = sum(rateBase(tms>0 & tms<responseWindow))
    % fine for now, but should calculate this on a trial to trial basis?
    % if time of response window is supposed to correspond to the duration of
    % first inhalation, which begins at a different time each trial

%these metrics are good checks but not as informative
nSpikesBase = sum(rateBaseHist(t>0 & t<responseWindow)*bin)/1000;
nSpikesBaseAvg = avgRateBase*responseWindow/1000;
    

rmax = 5*avgRateBase;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Now get stimulus-based properties     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rate = zeros(1,vp.nt);

    x    = zeros(1e5,1);    y    = zeros(1e5,1);
    nsp  = 0;
    kt   = 0;

    spkCounts_trials=[];
    for it = [stim(ks).in_tr]
        
        sp = sort(round(vertcat(tr(it).spikeTimes{un})));   %spike times for all requested units
        
        % If sType=odor, align to sniff; otherwise align to laser onset
        if stim(ks).odorConc >0  &&  strcmp(sType,'odor')
            %find first sniff after final valve open
            ii = find(tr(it).sniffZeroTimes(1,:)>tr(it).odorTimes(1)*1.009,1);
            if isempty(ii) || ii<1
                continue
            end
            t0 = tr(it).sniffZeroTimes(1,ii);  % t0 defined as inhalation onset
%             stim_str = sprintf('%s \n %1.1e', tr(it).odorName(1:min(length(tr(it).odorName),10)), tr(it).odorConc);   
        else
            t0 = tr(it).laserTimes(1);  % t0 defined as laser onset
%             stim_str = sprintf('Laser: \n%2.1f%s -%3.1fms', stim_intensity, stim_label, stim(ks).laserDur);
        end 
        kt = kt + 1;  % counting trials for current stim
        
        sp = sp((sp>t0+t1)&(sp<t0+t2))-t0; % spike times wrt t0 in window for current trial
        rate(sp-t1) = rate(sp-t1) + 1; % add 1 to position corresponding to ms when spikes occurred
        n  = numel(sp);
        
        % One way of counting the number of spikes every trial
        %  (how many spikes are from t0 to responseWindow)
        spkCounts_trials = [spkCounts_trials numel(sp((sp>0)&(sp<responseWindow)))];
        
        % Data for plotting raster
        x(nsp+(1:n)) = sp(:);
        y(nsp+(1:n)) = kt*ones(1,n);
        nsp = nsp + n;
        
        % Calculate ISI for each sp in trial
        ISI(kt).ISI_raw = nan(1,nt);
        if n > 1
            ISI(kt).ISI_raw(1,sp(2:end)-t1) = diff(sp);  % vector of ISIs
            ISI(kt).ISI_raw(1,sp(1)-t1) = inf;
            quickind = find(~isnan(ISI(kt).ISI_raw(1,:)));
            quicklook = ISI(kt).ISI_raw(1,quickind);   
        else
            quickind = [];
            quicklook = [];
        end  
    end
    spkCount = mean(spkCounts_trials)
    
    % Finished collecting spike data for each trial, so now trim raster and
    % reshape histogram vector
    x = x(1:nsp);
    y = y(1:nsp);
    
    rate = rate/kt;
    rateHist = sum(reshape(rate, bin, nt/bin),1)/(bin/1000); %units is Hz
    rateHist_smooth = smoothts(rateHist,'g',smooth.wsize,smooth.stdev);
    
        rate_hz = rate*1000;
    baseline=mean(rate_hz(tms<0));  
    pkThresh=baseline + 3*std(rate_hz(tms<0));
    [pkRate,pkInd]=findpeaks(rate_hz(tms>0),'MINPEAKHEIGHT',pkThresh,'NPEAKS',1);
    pkTime=tms(pkInd+find(tms>0,1))
    
    if isempty(pkTime)
        pkTime=nan;
        pkRate=nan;
    end

    % Another way of counting the number of spikes in response to stim
    nSpikes=sum(rate(tms>0 & tms<responseWindow));
    
    
%TODO: PLAY WITH DIFFERENT WAYS TO DEFINE THRESHOLD, SEE HOW MUCH LATENCY CHANGES
    resp_thr = 0.7 * pkRate;
    % resp_thr = pkThresh;
    isi_thr = 1/resp_thr*1000; % upper bound of ms since last spike for isi response
    
    % Pull stats on ISI from this response period
    ISItoplot = []; nsp = 0; firstISI = nan(kt,1); nSpkISI = zeros(kt,1);
    for jt = 1:kt
        isi_underthresh_t = find(ISI(jt).ISI_raw < isi_thr); %indices of spikes that followed the previous by less than threshold number of ms
        isi_underthresh_t = isi_underthresh_t(isi_underthresh_t>(0-t1) & isi_underthresh_t<(responseWindow-t1)); % just keep ones that are within response window
        
        if ~isempty(isi_underthresh_t)
            isi_underthresh = ISI(jt).ISI_raw(1,isi_underthresh_t);
            ISItoplot = [ISItoplot isi_underthresh];
            nSpkISI(jt) = numel(isi_underthresh_t); % note this does NOT count the spike that precedes the first under-threshold ISI
            
            firstISI(jt) = isi_underthresh_t(1) + t1;
        end
        
        % Save raster info for spikes used in ISI measures
        nn = numel(isi_underthresh_t); 
        xx(nsp+(1:nn)) = isi_underthresh_t +t1;   %%% histo of xx maybe related to pref sniff phase?? check with sniff stats, and with phase-triggered stimulus presentations
        yy(nsp+(1:nn)) = jt*ones(1,nn);
        nsp = nsp + nn;
    end

    latencyISI = round(mean(firstISI(isfinite(firstISI))))
    jitter = std(firstISI(isfinite(firstISI)));
    spkCountISI = mean(nSpkISI);
    avgISIResp = mean(ISItoplot);
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Save properties to response structure   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    vd.resp(ks).rateBase            = rateBase;
    vd.resp(ks).avgRateBase         = avgRateBase;
    vd.resp(ks).spkCountBase        = spkCountBase;
    vd.resp(ks).nSpikesBase         = nSpikesBase;
    vd.resp(ks).responseWindow      = responseWindow;
    %
    vd.resp(ks).rate                = rate;
    vd.resp(ks).spkCounts_eachTrial = spkCounts_trials;
    vd.resp(ks).spkCount            = spkCount;
    vd.resp(ks).spkCountErr         = std(spkCounts_trials);
    vd.resp(ks).nSpikes             = nSpikes;
    vd.resp(ks).maxFR               = pkRate;
    vd.resp(ks).latency             = pkTime;
    vd.resp(ks).latencyISI          = latencyISI;
    vd.resp(ks).jitterISI           = jitter;
    vd.resp(ks).spkCountISI         = spkCountISI;
    vd.resp(ks).avgISIResp          = avgISIResp;
    %
    vd.resp(ks).ntrial              = numel(stim(ks).in_tr);
    vd.resp(ks).stim                = stim(ks);
    vd.resp(ks).vp                  = vp;
    vd.resp(ks).raster.x            = x;
    vd.resp(ks).raster.y            = y;
    vd.resp(ks).inhLengths          = inhLengthsStim;
    vd.resp(ks).inhLengthsBase      = inhLenBase;
    vd.resp(ks).timeToFirstInhOdor  = timeToFirstInhOdor;
    %%%
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           Plot raster and PSTH            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ip = ip + 1;
    gf = figure(ifigRas);
    gs = subplot(2,numel(stim), ip);
    yt  = max([1, floor(kt/10)*10]);
    plot(x,y, '.', 'MarkerSize',7), hold on
%     	plot(xx,yy, '.g', 'MarkerSize',7), %only spikes with ISI below threshold
    plot([0,0], [0, kt+1], '--k'), hold off
    set(gs, 'XLim', [t1, t2], ...
        'YLim',      [0, kt+1], ...
        'YTick',     [0, yt], ...
        'XTick',     [-200,0,200], ...
        'FontSize',  10);
    title(stim(ks).stim_str, 'FontSize', 15, 'FontWeight', 'bold')
    
    
    gs = subplot(2, numel(stim), ip+numel(stim));
    plot(t, rateHist, ...
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
    
    %plot overlayed histograms if more than one stimulus
    if ks>1
        figure(ifigPsth)
        gs0 = subplot(1,1,1);
        plot(t, rateHist, ...
            'LineWidth', 2, ...
            'Color',     stim(ks).lineColor, ...
            'LineStyle', stim(ks).lineStyle),
        hold on
    end   
end

gf = figure(ifigRas);
sUnits=num2str(un(1));
for iu=2:numel(un)
    sUnits=sprintf('%s,%2d',sUnits,un(iu));
end
suptitle(sprintf('Mouse %s - session %3d - rec %s units %s',num2str(vp.mouse),str2num(vp.sess),vp.rec,sUnits));

%plot overlayed histograms if more than one stimulus
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Save the response data (the plot and the response struct) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(fn.fold_an_mouse, 'dir')
    mkdir(fn.fold_an_mouse)
end
if ~exist(fn.fold_an_sess, 'dir')
    mkdir(fn.fold_an_sess)
end

fileBase=[fn.basename_an sType '_units'];
unitstr = sprintf('%02d',un);
fileBase=[fileBase unitstr]

save(fullfile(fn.fold_an_sess,sprintf('%s_resp.mat',fileBase)), '-struct', 'vd','resp')
print(gf,'-depsc',fullfile(fullfile(fn.fold_an_sess,sprintf('%s_resp.eps',fileBase))))

resp=vd.resp;
end  %resp = view_rasters(sType,un)

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
stim = stim(in_sort);


for ks = 1:numel(stim)
    kt   = 0; stim(ks).lineColor = 'r';
    
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