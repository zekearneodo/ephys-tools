% Script for visualizing spikes aligned with the onset of stimuli
% uses data_managment_tools_032
% % 
% Updated by kp, 2014-05
%   can be used for recs with and without pulseOffset field
%   now contains isi analyses


function [vr,vd] = visualize_responses_kp_10(mouse,sess,rec,unit,sType)
global vp vd dm;  

%add the package of common files and folders to the path.
%here are the functions for file_naming, for instance, and any other
%functions that several script bundles will use.
%that is /baseFolder/ephysDataManagement/current/include
includePath=fullfile(fileparts(pwd),'current','include');
addpath(includePath);
% addpath('/Volumes/spikefolder')
% Functions of include that it uses:
%   - file_names(mouse,sess,rec,stat)

vr.visual_par_init     = @visual_par_init;
vr.make_stimuli_set    = @make_stimuli_set;
vr.view_stimuli_set    = @view_stimuli_set;
vr.analysis_file_names = @analysis_file_names;
% pp.ss_prep        = @ss_prep;
% pp.resampling     = @resampling;

if nargin < 5
    sType='laser'
end

dm = data_management_tools_032();
vp = visual_par_init(mouse,sess,rec,sType);
% vd = make_stimuli_set(vp);

view_stimuli_set(sType,unit);

end


function vp=visual_par_init(mouse,sess,rec,sType)
% visual_par_init(mouse,sess,rec,sType)
% sets the parameters for visualization, stored in global vp
% then calls function "vd = make_stimuli_set(vp)", which makes stimulus set
global dm vd

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

% define window to analyze and plot data
if strcmp(sType,'laser')
    vp.t1   =  -20;
    vp.t2   =  20;
    vp.respwin = 100;
else  % odor
    vp.t1   =  -200;
    vp.t2   =  400;
    vp.respwin = 200;
end
    
vp.nt  = vp.t2-vp.t1;
vp.bin = 1;
vp.t   = mean(reshape((vp.t1+1):vp.t2, vp.bin , vp.nt/vp.bin),1);
vp.smoothing = round(vp.nt/8); % bin size for gaussian smoothing of histo for plotting

vp.fn = dm.file_names(mouse, sess, rec);
q = load(vp.fn.trial);


%check if there is a trial correction function to be applied
trCFname=fullfile(vp.fn.fold_pr_sess,sprintf('%s_%s_trial_correct.m', vp.mouse,vp.sess));
if exist(trCFname)
    fprintf('Found trial correction function %s\n',trCFname);
    p = path;
%     p='/Volumes/spikefolder/';
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
if isfield(vp.tr,'pulseOnsetDelay')
    vp.tr_type = 3;
    vp.par = {'odorName', 'odorConc', 'laserDur', 'laserPower', 'pulseOnsetDelay'}
elseif isfield(vp.tr,'pulseOffset')        % just KPawakeM72_004
    vp.tr_type = 2;
    vp.par = {'odorName', 'odorConc', 'laserDur', 'laserPower', 'pulseGroup','pulseOffset'}
else                                       % old sessions
    vp.tr_type = 1;
    vp.par = {'odorName', 'odorConc', 'laserDur', 'laserAmp'}
end


%% select trials to look at, if you'd like

% % find indices of et tig trials
% for ii = 1:numel(vp.tr)
%     odorName{ii} = vp.tr(ii).odorName;
% end
% et = strfind(odorName,'ethyl tiglate');
% jj = 1;
% for ii = 1:numel(vp.tr)
%     if et{ii} == 1
%         et_in(jj) = ii;
%         jj = jj+1;
%     end
% end

% % find first trial after certain time in ms
% first_tr = find([vp.tr.start]>600000,1);

% select the trials you want to keep;
% in_tr=(first_tr:numel(vp.tr));
% vp.tr=vp.tr(in_tr)



vd = make_stimuli_set(vp);
end


function vd = make_stimuli_set(vp)
% ====================================================================================
% creates set of stimuli to plot, separating by properties, depending on 
% eg: odorstim,lightstim


for it=1:numel(vp.stimTypes)
    sType=vp.stimTypes{it}
    vd.(sType) = struct('list',[],'sort','')
end

switch vp.tr_type
    case 1
        stim = struct('odorName', '', 'odorConc', 0, 'laserDur', 0, 'laserAmp', 0, 'in_tr', [])
        tr=vp.tr;
    case 2
        stim = struct('odorName', '', 'odorConc', 0, 'laserDur', 0, 'laserPower', 0, 'in_tr', [], 'pulseOffset', [],'pulseGroup',[])
        tr=vp.tr;
    case 3
        stim = struct('odorName', '', 'odorConc', 0, 'laserDur', 0, 'laserPower', 0, 'in_tr', [], 'pulseOnsetDelay', [])
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
par = vp.par;
for it = 1:numel(tr)
    if tr(it).start == 0
        continue
    end
        
    tr(it).laserDur = diff(tr(it).laserTimes);
    
    new = 1; ks = 0;
    while new && (ks<ns)
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
        stim(ns).in_tr = it;
    end
end

stim = stim(2:end);
ns   = ns - 1;



% Remove from the list the stimuli for which there are less than 10 trials
keepStim=[];
for is=1:numel(stim)
    if numel(stim(is).in_tr)>10
        keepStim=[keepStim is];
    end
end
stim=stim(keepStim);


% Select property to group stimuli by
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
end 


function view_stimuli_set(sType,un)
global vd vp;

stim = vd.stim(vd.(sType).list);
sSort = vd.(sType).sort;

% Define order to present stimuli
if isnumeric(stim(1).(sSort))
    [~, in_sort] = sort(-[stim.(sSort)]);
else
    [~, in_sort] = sort({stim.(sSort)});
end

stim = stim(in_sort)


% Define properties for plotting stimuli
ns = numel(stim);
for ii = 1:ns
    stim(ii).lineStyle = '-';
    stim(ii).lineColor = 'k';
    if stim(ii).odorConc < 0.014
        stim(ii).lineStyle = '-';
        stim(ii).lineColor = 'r';
    end
end



%% =======================================================================
%   Plot raster and psth for defined units
%  =======================================================================

ifigRas=12;
% ifigPsth=2;
gf = figure(ifigRas); clf
ip = 0;
tr=vp.tr;
t1=vp.t1;
t2=vp.t2;
bin=vp.bin;
t=vp.t;
smoothing = vp.smoothing;
nt=vp.nt;
tms = t1+1:t2;
respwin = vp.respwin;
fn=vp.fn;
hold on


% Get baseline rate for odor and laser trials             
%    averaging sime size window around inhalation onset 
%    for the 3 blank sniffs before each trial start
rateBase=zeros(1,nt);
nb=0; inhLength_blank = []; inhLength_odor = [];
for kt = 1:numel(tr)
    trial=tr(kt);
    if isempty(trial.sniffZeroTimes) || isempty(trial.spikeTimes)
        continue
    end
    sp = sort(round(vertcat(trial.spikeTimes{un})));
    for nsb = 1:3       % now averaging the 3 sniffs before each trial
        iiBase = find(trial.sniffZeroTimes(1,:)>0,1) - nsb;
        if isempty(iiBase) || iiBase<1
            continue
        end
        t0Base = trial.sniffZeroTimes(1,iiBase);
        spBase = sp((sp>t0Base+t1)&(sp<=t0Base+t2))-t0Base;
        rateBase(spBase-t1)= rateBase(spBase-t1)+1;
        nb=nb+1;
        
        inhLength_blank = [inhLength_blank (trial.sniffZeroTimes(2,iiBase)-trial.sniffZeroTimes(1,iiBase))];
        if nsb==1 && trial.odorConc>0 
            inhLength_odor = [inhLength_odor (trial.sniffZeroTimes(2,iiBase+nsb)-trial.sniffZeroTimes(1,iiBase+nsb))];
        end
    end
end
figure; hold on
cdfplot(inhLength_blank)
cdfplot(inhLength_odor)
hold off

% spkBaseline = (1000/(0-t1))*sum(rateBase([t1:t2]<0))/nb;
foo=t1+1:t2; 
spkCountBase = sum(rateBase(foo>0 & foo<respwin))/nb;
rateBase = sum(reshape(rateBase, bin, nt/bin),1)/nb/(bin/1000);
avgRateBase = mean(rateBase);
nSpikesBase=sum(rateBase(t>0 & t<respwin)*bin)/1000;
rateBase = smoothts(rateBase,'g',smoothing);
    

%% Now, for each stimulus find spikes in analysis window for each trial
%  Store data for raster (x,y) and histogram (rate, reshaped)
for ks = 1:numel(stim)
%     1:numel(stim) [3 1 2]
    rate = zeros(1,nt);
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

    % Choose to plot only the first or last 20 trials, if you please
%     upperTrials=min(20,numel(stim(ks).in_tr));
%     stim(ks).in_tr=stim(ks).in_tr(end-upperTrials+1:end); 
%     stim(ks).in_tr=stim(ks).in_tr(1:upperTrials);
    
    
    spkCount=[];
    for it = [stim(ks).in_tr]
        % Get the spike times for all the requested units
        sp = sort(round(vertcat(tr(it).spikeTimes{un})));

        % If sType=odor, align to sniff; otherwise align to laser onset
        if stim(ks).odorConc >0 && strcmp(sType,'odor')
            ii = find(tr(it).sniffZeroTimes(1,:)>tr(it).odorTimes(1),1);  % first sniff after final valve open
            if isempty(ii) || ii<1
                continue
            end
            t0 = tr(it).sniffZeroTimes(1,ii);  % t0 defined as inhalation onset
            stim_str = sprintf('%s \n %1.1e', tr(it).odorName(1:min(length(tr(it).odorName),10)), tr(it).odorConc);
        
        else
            t0 = tr(it).laserTimes(1);  % t0 defined as laser onset
            stim_str = sprintf('Laser: \n%2.1f%s -%3.1fms', stim_intensity, stim_label, stim(ks).laserDur);
        end
        kt = kt + 1;  % count of trials for current stim
        
        sp = sp((sp>t0+t1)&(sp<t0+t2))-t0; % spike times wrt t0 in window for current trial
        rate(sp-t1) = rate(sp-t1) + 1; % add 1 to position corresponding to ms when spikes occurred
        n  = numel(sp);


        % One way of counting the number of spikes every trial
        %  (how many spikes are from t0 to t0+200)
        spkCount = [spkCount numel(sp((sp>0)&(sp<respwin)))]; % numel is number of trials
        
        % Data for plotting raster
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

    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx = sprintf('Stimulus #%i is %s',ks, stim_str)

    
    % Have collected data from spikes for each trial, so now trim raster and
    % reshape histogram vector
    x = x(1:nsp);
    y = y(1:nsp);

    spkCount2 = sum(rate(foo>0 & foo<respwin))/kt
    rate = sum(reshape(rate, bin, nt/bin),1)/kt/(bin/1000);
    rate_pl = smoothts(rate,'g',smoothing);
       
    
    % Calculate response measures (baseline, peak, etc)
    baseline = mean(rate(t<0));
    pkThresh = baseline + 3*std(rate(t<0));
    [pkRate,pkInd] = findpeaks([baseline rate(t>0)],'MINPEAKHEIGHT',pkThresh,'NPEAKS',1);
    pkTime = t(pkInd-1+find(t>0,1)-1);
%     spkBase = sum(rate(t<0))*(bin/1000);
    
    
    if isempty(pkTime)
        pkTime=0;
        pkRate=0;
        rmax = 40;
    else
        rmax = pkRate * 1.5;
    end
%     nSpikes=sum(rate(t>0 & t<respwin)*bin)/1000;
%     nSpikes_subtBase = nSpikes - nSpikesBase;
    rmax = 150;
    
    
    % Caluculate response period from PSTH
    resp_thr = 0.7 * pkRate; % might have to be higher for less strong responses
%     resp_thr = baseline + 3*std(rate(t<0));        

% PLAY WITH DIFFERENT WAYS TO DEFINE THRESHOLD, SEE HOW MUCH LATENCY CHANGES

%         resp_pd = [min(t(rate > resp_thr)) max(t(rate > resp_thr))];
    isi_thr = 1/resp_thr*1000 % upper bound for ms since last spike
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
    pkTime = pkTime
    latency = mean([ISI(isfinite([ISI(:).first])).first])
    jitter = std([ISI(isfinite([ISI(:).first])).first])
    avg_nSpks = mean([ISI(:).nSpks])
    avg_ISI = mean([ISI(isfinite([ISI(:).avgISI])).avgISI])
    
%     figure(3)
%     subplot(1,ns,ks)
%     hist(ISItoplot(ISItoplot<isi_thr),0:isi_thr/10:isi_thr)
%     hh = hist(ISItoplot(ISItoplot<isi_thr),0:isi_thr/10:isi_thr)
%     xlabel('ISI (ms), during response period only')
%     ylabel('Number of occurrences')
%     xlim([0 isi_thr])
%     ymax(2+max(hh))

    
    
    %% make response structure
    
    vd.resp(ks).avgRateBase   =   avgRateBase;
    vd.resp(ks).spkCount      =   spkCount;
%     vd.resp(ks).spkBase      =   spkBase;
    vd.resp(ks).spkBase       =   spkCountBase;
%     vd.resp(ks).spkBaseline  =   spkBaseline;
    vd.resp(ks).nSpikes       =   mean(spkCount);
    vd.resp(ks).nSpk_subtBase =   mean(spkCount)-spkCountBase;
    vd.resp(ks).nSpkErr       =   std(spkCount);
    vd.resp(ks).delay         =   pkTime;
    vd.resp(ks).maxFR         =   pkRate;
    vd.resp(ks).latencyISI    =   latency;
    vd.resp(ks).jitterISI     =   jitter;
    vd.resp(ks).nSpksISI      =   avg_nSpks;
    vd.resp(ks).respavgISI    =   avg_ISI;
    vd.resp(ks).stim          =   stim(ks);
    vd.resp(ks).ntrial        =   numel(stim(ks).in_tr);


    

    
    %%%
  
     %plots this psth
    ip = ip + 1;
    gf = figure(ifigRas);
    gs = subplot(2,numel(stim), ip);
    yt  = max([1, floor(kt/10)*10]); hold on
    plot(x,y, '.g', 'MarkerSize',7),
%         plot(xx,yy, '.', 'MarkerSize',7), %spikes with ISI below threshold
    plot([0,0], [0, kt+1], '--k'), hold off
    set(gs, 'XLim', [t1, t2], ...
        'YLim',      [0, kt+1], ...
        'YTick',     [0, yt], ...
        'XTick',     [t1:2:t2], ...
        'FontSize',  10);
    title(stim_str, 'FontSize', 10, 'FontWeight', 'bold')

    
    gs = subplot(2,numel(stim), ip+numel(stim));
    plot(t, rate_pl, 'g', ...
        'LineWidth', 4, ...
        'LineStyle', stim(ks).lineStyle),
%     'Color',     stim(ks).lineColor, ...
        hold on; box off
    %plot(t,ones(size(t))*baseline,'.k')
    %plot(t,ones(size(t))*pkThresh,'b')
    plot(t,rateBase,'LineStyle','--','LineWidth', 4,'Color',[.5,0.5,0.5])
%     plot(pkTime,pkRate,'k*')
    plot([0,0], [0, rmax], '--k'), hold off
    set(gs, 'XLim',     [t1, t2], ...
        'XTick',         [t1:200:t2], ...
        'YLim',          [0, rmax], ...
        'YTick',         0:100:rmax, ...
        'FontSize',      20);
    %         'XTick',         [t1:(t2-t1)/6:t2], ...
    
    %plots of psths altogether
%     if ks>1
%         figure(ifigPsth)
%         gs0 = subplot(1,1,1);
%         plot(t, rate, ...
%             'LineWidth', 2, ...
%             'Color',     stim(ks).lineColor, ...
%             'LineStyle', stim(ks).lineStyle),
%         hold on
%     end
    
 clear ISI xx yy   
end

    disp('oooooooooooooooooooooooooooooooooooooooooooooooo')


gf = figure(ifigRas);
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
   

% conci=[4 2 3]; affi = [3 1 5];
% figure; 
% plotyy(conci,[vd.resp(conci).latencyISI],conci,[vd.resp(conci).nSpksISI],'scatter')
% 
% figure; 
% plotyy(affi,[vd.resp(affi).latencyISI],[affi],[vd.resp(affi).nSpksISI],'scatter')

% figure; 
% plotyy([1:5],[vd.resp(:).latencyISI],[1:5],[vd.resp(:).nSpksISI])


% if ks>1
%     figure(ifigPsth)
%     gs0 = subplot(1,1,1);
%     subplot(gs0)
%     plot([0,0], [0, rmax], '--k'), hold off
%     set(gs0, 'XLim',     [t1, t2], ...
%         'XTick',         -200:100:200, ...
%         'YLim',          [0, rmax], ...
%         'YTick',         0:20:rmax, ...
%         'FontSize',      10);
% end
end
