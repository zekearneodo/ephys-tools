%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script for visualizing spikes aligned with the onset of stimuli

% outputs vr and response struct
% global vp contains the parameters for all the functions to call

%has functions:
% -visual_par_init
% -make_stimuli_set
% -view_rasters (default)
% -view_lfp (not called by default)

function [vr, resp] = visualize_responses_warped(mouse,sess,rec,unit,sType)
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

vp.warp = 1;  % 1 if want to plot in warped space, 0 if not
vp.sort = 1;  % 1 if want to sort trials by length of 1st inhalation, 0 if order presented


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
    vp.warp = 0;
    vp.sort = 0;
else
    vp.t1  = -400;
    vp.t2  =  1500;
    vp.responseWindowDefault=400;
end

% vp.nt  = vp.t2-vp.t1;
vp.bin = 20;
% vp.t   = mean(reshape((vp.t1+1):vp.t2, vp.bin , vp.nt/vp.bin),1);
% vp.smooth.wsize = round(vp.nt/8);   %window size for gaussian smoothing of histo for plotting
% vp.smooth.cutoff = 100;   %cutoff for the gaussian smoothing
% vp.smooth.stdev = q.info.rec(nrec).sampling_freq/(2*pi*vp.smooth.cutoff); %std for gaussian smoothing of histo for plotting

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
clear q

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
            same_stim = same_stim & (tr(it).(par{ip}) == stim(ks).(par{ip}));
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

try
    odorInfo(1).name      = voPar.odor;
    try
        odorInfo(1).flows     = [voPar.AirFlow_1 ;  voPar.NitrogenFlow_1];
    catch
        odorInfo(1).flows     = [voPar.air ;  voPar.nitrogen];
    end
    
    if isfield(voPar,'dillution')
        dil=voPar.dillution;
    else
        dil=1;
    end
    odorInfo(1).dillution = [dil];
    odorInfo(1).vialId    = [voPar.vial];
    odorInfo(1).vialConc  = [voPar.vialconc];
catch
    warning('There were errors trying to get odor info, it might be incomplete');
end

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

    tr  = vp.tr;
    bin = vp.bin;
    t1  = vp.t1;
    t2 = vp.t2;
%     if strcmp(sType,'odor') && vp.warp
%         FVdur = mode(diff([tr(diff([tr.odorTimes])>0).odorTimes]));
%         t2=max(max([tr.odorTimes]));  %% end of t window set to max of FV off
%         clipt=rem(t2,bin);
%         t2=t2-clipt;
%     else
%         t2 = vp.t2;
%     end
    if rem(t2,bin) ~=0
        t2=t2-rem(t2,bin);
    end
    nt  = t2-t1;
    tms = t1+1:t2;  %time from t1:t2 in ms, not reshaped
    t   = mean(reshape(tms, bin , nt/bin),1);  %time from t:t2, reshaped according to bin size
    vp.smooth.wsize = round(nt/8);   %window size for gaussian smoothing of histo for plotting
    vp.smooth.cutoff = 100;   %cutoff for the gaussian smoothing
    vp.smooth.stdev = 19358/(2*pi*vp.smooth.cutoff); %std for gaussian smoothing of histo for plotting
    smooth = vp.smooth;
    fn  = vp.fn;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%            Get baseline rate              %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Averaging spikes in same size window (t1:t2) around inhalation onsets
    % for a few blank sniffs before each trial start.
    bad_trials=[];
    rateBase=zeros(1,nt);   spBase_byTrial=zeros(1,nt);
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
        %%% Can use this as a quick way to cut out trials where nonstationary
        %%% units were lost
%         if numel(sp>-500 & sp<1500)<3
%             bad_trials = [bad_trials kt];
%             continue
%         end
        nWindow=0;
        oneMoreWindow=1;
        
        while oneMoreWindow
            % This window is the first that that is part of the stimulus period
            if nWindow==0
                if trial.odorConc>0
                    iFirstSniff = find(trial.sniffZeroTimes(1,:)>trial.odorTimes(1),1,'first');
                    %                 iFirstSniff = iFirstSniff+1;
                    t0Base = trial.sniffZeroTimes(1,iFirstSniff);  %time of first sniff after stimulus
                    timeToFirstInhOdor = [timeToFirstInhOdor t0Base-trial.odorTimes(1)];
                    
                elseif trial.laserAmp>0 && any(strcmpi(trial.odorName,{'none','empty','dummy'}))
                    iFirstSniff = find(trial.sniffZeroTimes(1,:)<trial.laserTimes(1)*1.009,1,'last'); %t0 is TRIAL start (not FV or laser), could correct based on stim type
                    t0Base = trial.sniffZeroTimes(1,iFirstSniff);  %time of first sniff after stimulus
                    
                    laserTimeAfterInhDetected = [laserTimeAfterInhDetected round(1.009*trial.laserTimes(1)-t0Base)];
                else
                    break
                end

                if isempty(t0Base)
                    break
                end
                %%%%  Respiration statistics:
                % get duration of first inh after stim onset
                inhLenStim = [inhLenStim diff(trial.sniffZeroTimes(:,iFirstSniff))];
                % get durations of the three sniff cycles and inhalations
                %   after stimulus onset
                for nS=nSniffs-(0:nSniffs-1)
                    try
                        sniffLengths(kt,1:nS) = diff(trial.sniffZeroTimes(1,iFirstSniff:iFirstSniff+nS));
                        inhLengthsStim(kt,1:nS+1) = diff(trial.sniffZeroTimes(:,iFirstSniff:iFirstSniff+nS));
                        break
                    catch
                        continue
                    end
                end
                % get all pre-stimulus inhale durations
                inhLenBase = [inhLenBase diff(trial.sniffZeroTimes(:,1:iFirstSniff-1))];
                
            % Now add up the spikes in as many non-overlapping windows as fit in the pre-stimulus period
            else
                iFirstSniff = find(trial.sniffZeroTimes(1,:)<t0Base-nt,1,'last');
                
                if isempty(iFirstSniff) || iFirstSniff<1
                    break
                elseif trial.sniffZeroTimes(2,iFirstSniff)>0
                    break
                end

                nb=nb+1; %keep track of the number of *segments* added to baseline

                t0Base = trial.sniffZeroTimes(1,iFirstSniff);
                spBase = sp((sp>t0Base+t1)&(sp<=t0Base+t2))-t0Base;
                spBase_byTrial(nb,spBase-t1) = 1;   
%                 spBase_tr_10ms(nb,:) = sum(reshape(spBase_trials(nb,:),10,nt/10));
                rateBase(spBase-t1) = rateBase(spBase-t1)+1;
                
                %%% Baseline for cumulative count latency measure
                spBase_cC = spBase;
                %  spBase_cC = spBase(spBase>0);
                for cc = 1:numel(spBase_cC)
                    cumulativeCountBase(nb,(spBase_cC(cc)-t1):nt) = cc;
                end
                
                if t0Base-nt < min(trial.sniffZeroTimes(1,:))
                    oneMoreWindow=0;
                end

            end
            nWindow=nWindow+1;
        end
    end
    
    rateBase=rateBase/nb; %divide spk times collected by number of segments used
    
    
% %  %%% Pieces to look at variance of spiking in future
% %     base_std = std(spBase_byTrial);
% %     base_var = var(spBase_byTrial);
% %     var_ratio = base_var./rateBase;
% %     avg_var_ratio = mean(var_ratio(isfinite(var_ratio)));
% %     
% %     threshRate_ms = (rateBase*1000 + 2000*base_std);
% % 
% %     figure; hold on
% %     plot(tms,rateBase*1000,'-','LineWidth',3)
% %     plot(tms,threshRate_ms,'r--')
% %     plot(tms,base_std*1000,'--','Color',[0.5 0.5 0.5])

    
    %% This shows if we had an issue with laser triggering
    laseroffset=laserTimeAfterInhDetected;
    %%
    
% % %     %%%%% Plot respiration statistics
% % %     figure(1); hold on
% % %     subplot(3,1,1)
% % %     cdfplot(inhLenBase)
% % %     xlim([0 300]);
% % %     title('no odor inhale duration')
% % %     subplot(3,1,2)
% % %     cdfplot(inhLenStim)
% % %     xlim([0 300]);
% % %     title('first inhale with stim duration')
% % %     subplot(3,1,3)
% % %     cdfplot(timeToFirstInhOdor)
% % %     xlim([0 1000]);
% % %     title('latency to first inh after FV on')
% % %     hold off
% % %       
    
    
    %% =======================================================================
    %%  For each stimulus find spikes in analysis window for each trial
    %   Store data for raster (x,y) and histogram (rate, reshaped)
    
    % Set up the figure
    ifigRas=figure;
    scrsz = get(0,'ScreenSize');
    set(ifigRas,'Position',[1 (scrsz(4)/2) 3*scrsz(3)/4 scrsz(4)/2],...
        'Nextplot','add');
    nSubPlots=double(2*numel(stim)); hS=zeros(1,nSubPlots);
    figRows=2;  figCols=ceil(nSubPlots/figRows);
    for iSubFig=1:nSubPlots
        hS(iSubFig)=subplot(figRows,figCols,iSubFig);
        set(hS(iSubFig),'Nextplot','add');
    end
    
    ip  = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%    First, calculate baseline properties   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Crop extra 1ms bins so time and rate vectors can be reshaped
%     extrams = rem(nt,bin);
%     nt  =  t2-t1-extrams;
%     tms =  t1:t2-1-extrams;  %time from t1:t2 in ms, not reshaped
%     t   =  mean(reshape(tms,bin ,nt/bin),1);
    
    %         rateBase_warped = rateBase_warped(1:end-extrams);
    rateBaseHist = sum(reshape(rateBase, bin, nt/bin),1)/(bin/1000); %units is Hz
    rateBaseHist_smooth = smoothts(rateBaseHist,'g',smooth.wsize,smooth.stdev);
    avgRateBase  = mean(rateBase*1000) %units is Hz
    stdRateBase  = std(rateBase*1000);
    rateBase_std_bin = std(reshape(rateBase*1000, bin, nt/bin),1);
    
    
    
    for ks = 1:numel(stim)
        
        % Build subplot titles with relevant stimulus info
        switch 1
            case (stim(ks).odorConc>0) && ~(stim(ks).laserDur>0)
                stim_str = sprintf('%s \n %1.1e', stim(ks).odorName(1:min(length(stim(ks).odorName),10)), stim(ks).odorConc);
            case ~(stim(ks).odorConc>0) && (stim(ks).laserDur>0)
                stim_str = sprintf('Laser: \n%2.1fmV -%3.1fms', stim(ks).laserAmp, stim(ks).laserDur);
            case (stim(ks).odorConc>0) && (stim(ks).laserDur>0)
                stim_str = sprintf('Odor+Laser: %s \n %2.1fmV -%3.1fms',stim(ks).odorName(1:min(length(stim(ks).odorName),10)),stim(ks).laserAmp, stim(ks).laserDur);
        end
        sprintf('=======================================\nStim #%i \n%s',ks, stim_str)
        

        
        % Define responseWindow (time after stim you count n spikes)
        if ~any(strcmpi(stim(ks).odorName,{'none','dummy'}))
            % vp.responseWindow=mean(inhLenStim);
            vp.responseWindow=vp.responseWindowDefault;
        else
            vp.responseWindow=vp.responseWindowDefault;
        end
        
        responseWindow = vp.responseWindow;
        spkCountBase   = sum(rateBase(tms>0 & tms<responseWindow));
        nSpikesBaseAvg = avgRateBase*responseWindow/1000;
        
        
        
        
        rmax = 40;

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%     Now get stimulus-based properties     %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        rate = zeros(1,nt);
        
        FULLpreStim_resp_times=nan(2e2,30);
        FULLpostStim_resp_times=nan(2e2,50);

        %         cumulativeCount = zeros(1e3,nt);
        x    = zeros(1e5,1);    y    = zeros(1e5,1);
        nsp  = 0;
        kt   = 0;   wkt = [];    kt_skip=[];
        spkCount=[];	
        
        for it = [stim(ks).in_tr]
            if isempty(tr(it).start)
                continue
            end
            %%For nonstationary units, can skip trials with no spike events
            %   (bad_trials defined while getting baseline rate)
%             if any(bad_trials==it)
%                 aaa=573;
%                 continue
%             end

            % Concatenate spike times from all requested units
            sp = sort(round(vertcat(tr(it).spikeTimes{un})))';
            
            % Odor: reference spike times to first inhalation
            % else: align to onset of laser
            if stim(ks).odorConc >0 && strcmp(sType,'odor')
                % ii is index of first inhalation after final valve open
                ii = find(tr(it).sniffZeroTimes(1,:)>tr(it).odorTimes(1)*1.009);
                
                if isempty(ii) || numel(ii)<2
                    continue
                end
                ii=ii(1);
                t0 = tr(it).sniffZeroTimes(1,ii); %time reference for all events of trial
                
                if t0>tr(it).odorTimes(2)*1.009
                    continue
                end
                
                
                
                kt = kt + 1;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Build matrix of respiration times, referenced to t0
                keepthese=[];
                prephases=find(tr(it).sniffZeroTimes(1,:)>t1,1):ii-1;
                postphases=ii:find(tr(it).sniffZeroTimes(1,:)>t2,1);
                if numel(postphases)<3 || numel(prephases)<1
                    warning('May not be able to analyze this stimulus with warping')
                    kt_skip = [kt_skip kt];
                    continue
                else
                wkt=[wkt kt];
                for iph=1:numel(prephases)
                    keepthese(1,2*iph-1:2*iph)=tr(it).sniffZeroTimes(:,prephases(iph))'-t0;
                end
                FULLpreStim_resp_times(numel(wkt),(end+1-numel(keepthese)):end)=keepthese;

                for iph=1:numel(postphases)-1
                    FULLpostStim_resp_times(numel(wkt),2*iph-1:2*iph)=tr(it).sniffZeroTimes(:,postphases(iph))'-t0;
                end
                end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            elseif strcmp(sType,'laser')
                t0 = round(tr(it).laserTimes(1)*1.009);
                kt = kt + 1;
            end
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Collect spikes in requested window for raster and psth
            sp = sp((sp>t0+t1)&(sp<t0+t2))-t0;
            rate(sp-t1) = rate(sp-t1) + 1;
            n  = numel(sp);
            
            x(nsp+(1:n)) = sp(:);
            y(nsp+(1:n)) = kt*ones(1,n);
            nsp = nsp + n;
            
            %Keeping track of nSpikes in defined response window
            spkCount=[spkCount numel(sp((sp>0)&(sp<responseWindow)))];
            
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
            
            % Cumulative spike count (for another measure of latency)
            sp_cC = sp;
            %         sp_cC = sp(sp>0);
            n_cC = numel(sp_cC);
            for cc = 1:n_cC
                cumulativeCount(kt,(sp_cC(cc)-t1):nt) = cc;
            end
        end
        
        
        spkCount = mean(spkCount);
        
        % Finished collecting spike data for each trial, so now trim raster and
        % reshape histogram vector
        x = x(1:nsp);
        y = y(1:nsp);
%         cumulativeCount = cumulativeCount(1:kt,:);
        
        rate = rate/kt;
        rate_hz = rate*1000;
        
        rateHist = sum(reshape(rate, bin, nt/bin),1)/(bin/1000); %units is Hz
%         rateHist_smooth = smoothts(rateHist,'g',smooth.wsize,smooth.stdev);

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Normalize spike times to respiration phase lengths
        
            % Find min number columns (sniffs) with no nans
                %%%Could change to allowing one or two nans (and leaving
                %%%out those trials) so have more sniffs to analyze
            nSniff_pre = numel(find(sum(isnan(FULLpreStim_resp_times(1:numel(wkt),:)),1)==0));
            nSniff_post = numel(find(sum(isnan(FULLpostStim_resp_times(1:numel(wkt),:)),1)==0));
            resp_t0 = nSniff_pre+1;
            
            FULLpreStim_resp_times =FULLpreStim_resp_times(1:numel(wkt),end+1-nSniff_pre:end);
            FULLpostStim_resp_times =FULLpostStim_resp_times(1:numel(wkt),1:nSniff_post);
            
            respirationTimes = [FULLpreStim_resp_times FULLpostStim_resp_times];
            norm_respirationTimes = round(mean(respirationTimes,1));
            max_respCycles = 3;
            xwin = [norm_respirationTimes(resp_t0-min(nSniff_pre,2)) norm_respirationTimes(resp_t0+min(nSniff_post-1,max_respCycles*2+1))];
            tw = xwin(1):xwin(2);
        
        if vp.warp && strcmp(sType,'odor')
            badTrials=[];
            x_warped=[];
            y_warped=[];
            rate_warped=zeros(1,xwin(2)-xwin(1)+1);
            for itr=wkt 
                iw=find(wkt==itr);
                spikeTimes_thisTrial = x(y==itr);
                normalizedSpikes=[];
                
                for iphase = find(norm_respirationTimes<xwin(2) & norm_respirationTimes>=xwin(1))
                phaseBounds=respirationTimes(iw,iphase:iphase+1);
                phaseLength=diff(phaseBounds);
                ind_toWarp = find(spikeTimes_thisTrial>phaseBounds(1) & spikeTimes_thisTrial<phaseBounds(2));
                
                if ~isempty(ind_toWarp)
                    % Check for possible missed sniffs and remove those trials
%                     if phaseLength>1500
%                         badTrials=[badTrials find(phaseLength>1000)]
%                         continue
%                         %phaseLengths(badTrials) = nan;
%                     end
                    
                    % Set up warping
                    norm_phaseLength(iphase) = norm_respirationTimes(iphase+1)-norm_respirationTimes(iphase);
                    normalizedTime = linspace(norm_respirationTimes(iphase),norm_respirationTimes(iphase+1),norm_phaseLength(iphase)+1);
                    
                    spikeTimes_toWarp=spikeTimes_thisTrial(ind_toWarp);
%                       warped spk time within phase = (sp-t1)*(tm/(t2-t1))
                    normalizedSpikes = [normalizedSpikes; round((spikeTimes_toWarp-phaseBounds(1)) .* (norm_phaseLength(iphase) / phaseLength)) + norm_respirationTimes(iphase)];

                    % Normalize the time of each spike to the resp
                    %  phase length of this trial
%                     for itr_x = 1:numel(spikeTimes_toWarp)
%                         norm_t_thisSniff = linspace(phaseBounds(1),phaseLength+phaseBounds(1),norm_phaseLength(iphase)+1);
%                         [~,closest_normTime] = min(abs(norm_t_thisSniff-spikeTimes_toWarp(itr_x)));
%                         normalizedSpikes1 = [normalizedSpikes1; normalizedTime(closest_normTime)];
%                         try
%                             rate_warped1(1,normalizedTime(closest_normTime)-xwin(1)) = 1+ rate_warped1(1,normalizedTime(closest_normTime)-xwin(1));
%                         catch
%                             continue
%                         end
%                         
%                     end
                end
                end
                rate_warped(normalizedSpikes-xwin(1)) = 1 + rate_warped(normalizedSpikes-xwin(1));

                x_warped = [x_warped; normalizedSpikes];
                y_warped = [y_warped; iw.*ones(length(normalizedSpikes),1)];
            end


            
            %%% Create histogram of warped data
            extrams = rem(length(rate_warped),bin);
            tw=tw(1:end-extrams);
            tw = sum(reshape(tw, bin, length(tw)/bin),1)/bin;
            rate_warped = rate_warped(1:end-extrams)/numel(wkt);
            rateHist_warped = sum(reshape(rate_warped, bin, length(rate_warped)/bin),1)/(bin/1000);
            
            %%% !!! This method of finding latency (single static threshold)
            %%%     does not work for phasic cells  !!!
            warp_thresh = mean(rateHist_warped(tw<0)) + 3*std(rateHist_warped(tw<0));
%             warp_thresh(warp_thresh<1) = 1;
            [firstPk_warp_rate,pkInd_w]= findpeaks(rateHist_warped(tw>0),'NPeaks',1,'MinPeakHeight',warp_thresh);
            firstPk_warp_time = tw(tw>0);
            firstPk_warp_time = firstPk_warp_time(pkInd_w);
            
%             figure;
%             subplot(2,1,1)
%             plot(x,y,'b.')
%             xlim([t1 t2])
%             subplot(2,1,2)
%             plot(x_warped,y_warped,'g.')
%             hold on
%             plot(tw,rateHist_warped)
%             xlim([norm_respirationTimes(1) norm_respirationTimes(end)])
% % %             fill(norm_respirationTimes)
% %            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Create vectors for plotting inhalation periods in grey
        
        %%% COULD USE TRY/CATCH HERE TO EXPEDITE?
        
        [inh1_length,inh1_idx]=sort(respirationTimes(:,resp_t0+1));
        if ~vp.sort && ~vp.warp
            inh_plot_x = [0 0 respirationTimes(:,resp_t0+1)'];
            inh_plot_y = [length(inh1_idx) 0 1:length(inh1_idx)];
            if nSniff_post>3
            inh2_plot_x = [respirationTimes(1,resp_t0+2) respirationTimes(:,resp_t0+2)' flipud(respirationTimes(:,resp_t0+3))' respirationTimes(1,resp_t0+3)];
            inh2_plot_y = [0:size(respirationTimes,1) size(respirationTimes,1):-1:0];
            if nSniff_post>5
            inh3_plot_x = [respirationTimes(1,resp_t0+4) respirationTimes(:,resp_t0+4)' flipud(respirationTimes(:,resp_t0+5))' respirationTimes(1,resp_t0+5)];
            inh3_plot_y = [0:size(respirationTimes,1) size(respirationTimes,1):-1:0];
            end
            end
        elseif vp.sort && ~vp.warp
            inh_plot_x = [0 0 inh1_length'];
            inh_plot_y = [length(inh1_idx) 0 1:length(inh1_idx)];
            if nSniff_post>3
            inh2_plot_x = [respirationTimes(1,resp_t0+2) respirationTimes(:,resp_t0+2)' flipud(respirationTimes(:,resp_t0+3))' respirationTimes(1,resp_t0+3)];
            inh2_plot_y = [0:size(respirationTimes,1) size(respirationTimes,1):-1:0];
            if nSniff_post>5
            inh3_plot_x = [respirationTimes(1,resp_t0+4) respirationTimes(:,resp_t0+4)' flipud(respirationTimes(:,resp_t0+5))' respirationTimes(1,resp_t0+5)];
            inh3_plot_y = [0:size(respirationTimes,1) size(respirationTimes,1):-1:0];
            end
            end
        elseif vp.sort && vp.warp
            inh_plot_x = [0 0 norm_respirationTimes(resp_t0+1).*ones(1,length(inh1_idx)+1)];
            inh_plot_y = [length(inh1_idx) 0 0:length(inh1_idx)];
            if nSniff_post>3
            inh2_plot_x = [norm_respirationTimes(resp_t0+2).*ones(1,length(inh1_idx)+1) norm_respirationTimes(resp_t0+3).*ones(1,length(inh1_idx)+1)];
            inh2_plot_y = [0:size(respirationTimes,1) size(respirationTimes,1):-1:0];
            if nSniff_post>5
            inh3_plot_x = [norm_respirationTimes(resp_t0+4).*ones(1,length(inh1_idx)+1) norm_respirationTimes(resp_t0+5).*ones(1,length(inh1_idx)+1)];
            inh3_plot_y = [0:size(respirationTimes,1) size(respirationTimes,1):-1:0];
            end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Sort trials in raster by length of inhalation
        
        if vp.sort
            y_sorted=[];    x_sorted=[];
            for iin=1:length(inh1_idx)
                this_inh=inh1_idx(iin);
                if vp.warp
                    this_y = find(y_warped==this_inh);
                    y_sorted(end+1:end+length(this_y)) = iin;
                    x_sorted(end+1:end+length(this_y)) = x_warped(this_y);
                else
                    this_y = find(y==this_inh);
                    y_sorted(end+1:end+length(this_y)) = iin;
                    x_sorted(end+1:end+length(this_y)) = x(this_y);
                end
            end
        end
     
%         ccThresh = mean(cumulativeCountBase,1) + 2.05*std(cumulativeCountBase,1);
        %%%%%%%%%
        %%% Plot cumulative count of spikes
        %     figure; hold on
        % %     plot(tms,cumulativeCountBase(1:5:250,:),'k')
        %     plot(tms,cumulativeCount,'r')
        %     plot(tms, mean(cumulativeCount,1), 'k', 'LineWidth', 4)
        %     plot(tms, mean(cumulativeCountBase,1), 'Color', [0.5 0.5 0.5], 'LineWidth', 4)
        %     plot(tms, ccThresh, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 2)
        %     title(['Cumulative spike count -- ',stim(ks).stim_str])
        
%         combinedCount = vertcat((mean(cumulativeCount,1)),ccThresh);
%         diff_combinedCount = diff(combinedCount,1);
%         diff_combinedCount(diff_combinedCount>0) = 1;
%         diff_combinedCount(diff_combinedCount<0) = 0;
%         latency_cumulativeCount = find(diff(diff_combinedCount),1) + t1;
        %%%%%%%%
        
                
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Another latency calculation:
        %%%     Interpolate point where rateHist crosses time-varying threshold
        
        ind_crossing = [];  latency_interpolated = [];  rate_linfit = [];
        
        %%% i) identify first 10ms bin after t0 where rate>rateBase
        rateThresh_bin = rateBaseHist + 3*rateBase_std_bin;
        
        overThresh_ind = length(t(t<0))+find(rateHist(t>0)>rateThresh_bin(t>0),1);
        
        %             figure; hold on
        %             plot(t,rateBase_bin,'Color',[0.5 0.5 0.5])
        %             plot(t,rateThresh_bin,'c')
        %             plot(t,rateHist,'r')
        
        %%% ii) interpolate to find where rate exceeds threshold
        if ~isempty(overThresh_ind)
            thresh_linfit = linspace(rateThresh_bin(overThresh_ind-1),rateThresh_bin(overThresh_ind),1000);
            rate_linfit = linspace(rateHist(overThresh_ind-1),rateHist(overThresh_ind),1000);
            
            ind_crossing = find(diff(sign(thresh_linfit-rate_linfit)));
            latency_interpolated = round(t(overThresh_ind-1) + bin/1000*ind_crossing)
        else
            overThresh_ind = nan;
            latency_interpolated = nan
            rate_linfit = nan;
            ind_crossing = 1;
        end
        %             plot(time_crossing,rate_linfit(ind_crossing),'*g')
        
        %%% iii) now find first significant suppression
        rateThresh_bin = rateBaseHist - 3*rateBase_std_bin;
        underThresh_ind = length(t(t<0))+find(rateHist(t>0)<rateThresh_bin(t>0),1);
        
        if ~isempty(underThresh_ind)
            thresh_linfit_inh = linspace(rateThresh_bin(underThresh_ind-1),rateThresh_bin(underThresh_ind),1000);
            rate_linfit_inh = linspace(rateHist(underThresh_ind-1),rateHist(underThresh_ind),1000);
            
            ind_inh = find(diff(sign(thresh_linfit_inh-rate_linfit_inh)));
            inhibition_interpolated = round(t(underThresh_ind-1) + bin/1000*ind_inh)
        else
            underThresh_ind = nan;
            inhibition_interpolated = nan
            rate_linfit_inh = nan;
            ind_inh = 1;
        end
        
        %%
        %%%%% Next:
        %%%%%   - ) Repeat interpolation for suppression
        %%%%%   ii) Trial by trial
        %%%%%   - ) Vary parameters (n*std, bin size, phasic vs flat cells)
        %%%%%   - ) Find fist and max peaks
        %%%%%   - ) Find end of and width of excitation
        %%%%%   - ) Sort trials by inhalation length
        %%%%%  vii) Keep only trials in the middle 80% of inhalation lengths
        %%%%%
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%  Find first peak and max peak of the binned histogram
        
        % Manual find peaks code, so that threshold can vary with time
        t_stim=t(t>0);                  rateBase_stim = rateBaseHist(t>0);
        rateHist_stim=rateHist(t>0);    rateBase_std_stim=rateBase_std_bin(t>0);
        
        manPks_in = 1+find(diff(sign(diff(rateHist_stim)))==-2);
        manPks = rateHist_stim(manPks_in) > (3*rateBase_std_stim(manPks_in) +rateBase_stim(manPks_in));
        manPks_in = manPks_in(manPks);
        
        manPks_time = t_stim(manPks_in);
        manPks_rate = rateHist_stim(manPks_in);
        
        if isempty(manPks_rate)
            firstPk_time=nan
            firstPk_rate=nan;
            maxPk_time=nan
            maxPk_rate=nan;
        else
            firstPk_time = manPks_time(1)
            firstPk_rate = manPks_rate(1);
            maxPk_time = t_stim(find(rateHist_stim==max(manPks_rate),1))
            maxPk_rate = max(manPks_rate);
        end
        
        
        
        % Another way of counting the number of spikes in response to stim
        nSpikes=sum(rate(tms>0 & tms<responseWindow));
        
        
        %TODO: PLAY WITH DIFFERENT WAYS TO DEFINE THRESHOLD, SEE HOW MUCH LATENCY CHANGES
        %             resp_thr = 0.5 * maxPk_rate;
        resp_thr = 3*mean(rateBase_std_stim) + mean(rateBase_stim);
        %     resp_thr = 3*stdRateBase + avgRateBase;
        isi_thr = 1/resp_thr*1000; % upper bound of ms since last spike for isi response
        
        % Pull stats on ISI from this response period
        ISItoplot = []; nsp = 0; firstISI = nan(kt,1); nSpkISI = zeros(kt,1);
        for jt = 1:numel(ISI)
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
        if strcmp(sType,'odor')
            %                 firstISI_sorted = firstISI(in_inh1);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%   Save properties to response structure   %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        vd.resp(ks).rateBase            = rateBase;
        vd.resp(ks).avgRateBase         = avgRateBase;
        vd.resp(ks).spkCountBase        = spkCountBase;
        vd.resp(ks).responseWindow      = responseWindow;
        %
        vd.resp(ks).rate                = rate;
        vd.resp(ks).spkCounts_eachTrial = spkCount;
        vd.resp(ks).spkCount            = spkCount;
        vd.resp(ks).spkCountErr         = std(spkCount);
        vd.resp(ks).nSpikes             = nSpikes;
        vd.resp(ks).rateMax             = maxPk_rate;
        vd.resp(ks).latencyMax          = maxPk_time;
        vd.resp(ks).rateFirstPk         = firstPk_rate;
        vd.resp(ks).latencyFirstPk      = firstPk_time;
        vd.resp(ks).latencyISI          = latencyISI;
        vd.resp(ks).jitterISI           = jitter;
        vd.resp(ks).spkCountISI         = spkCountISI;
        vd.resp(ks).avgISIResp          = avgISIResp;
        vd.resp(ks).latency             = latency_interpolated;
        vd.resp(ks).inhibitionDetected  = inhibition_interpolated;
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
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Plot raster in upper subplot
        ip = ip + 1;
        gf = figure(ifigRas);
        gs = hS(ip);
        subplot(gs); hold on
        if strcmp(sType,'odor')
        fill(inh_plot_x,inh_plot_y,[0.95 0.95 0.95],'LineStyle','-','LineWidth',0.01)
        if nSniff_post>3
        fill(inh2_plot_x,inh2_plot_y,[0.95 0.95 0.95],'LineStyle','-','LineWidth',0.01)
        if nSniff_post>5
        fill(inh3_plot_x,inh3_plot_y,[0.95 0.95 0.95],'LineStyle','-','LineWidth',0.01)
        end
        end
        end
        
        plot(x_sorted,y_sorted, 'b.', 'MarkerSize',6)
        set(gs, 'XLim', xwin,...
            'YLim',      [0 kt+1],...
            'YTick',     [0 numel(wkt)],...
            'XTick',     [fliplr(0:-200:xwin(1)) 200:200:xwin(2)],...
            'FontSize',  10);
        %plot(firstISI_sorted,1:length(firstISI_sorted), 'g.', 'MarkerSize',6)
        plot([0,0], [0, kt+1], '--k'), hold off
        title(gs,stim_str,'FontSize', 14)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot histogram in lower subplot
        gs = hS(ip+numel(stim));
        subplot(gs); hold on
        if ~vp.warp
        plot(t, rateHist, ...
            'LineWidth', 1, ...
            'Color',     stim(ks).lineColor, ...
            'LineStyle', stim(ks).lineStyle)
        hold on
        set(gs, 'XLim',     [t1, t2], ...
            'XTick',         [t1:500:0 500:500:t2], ...
            'YLim',          [0, rmax], ...
            'YTick',         0:50:rmax, ...
            'FontSize',      12);
        plot(maxPk_time,maxPk_rate,'k*')
        plot(firstPk_time,firstPk_rate,'k*')
        % plot(latencyISI,resp_thr,'c*','MarkerSize',6)
        plot(latency_interpolated,rate_linfit(ind_crossing),'*b','MarkerSize',8)
        plot(inhibition_interpolated,rate_linfit_inh(ind_inh),'ob','MarkerSize',6)
        elseif vp.warp
            plot(tw,rateHist_warped, ...
                'LineWidth', 1, ...
                'Color',     stim(ks).lineColor, ...
                'LineStyle', stim(ks).lineStyle)
            plot(firstPk_warp_time,firstPk_warp_rate,'k*')
%             plot(firstPk_time,firstPk_rate,'b*')
            set(gs,'XLim',   xwin,...
                'XTick',     [xwin(1):500:0 500:500:xwin(2)],...
                'YLim',      [0 rmax],...
                'YTick',     0:50:rmax,...
                'FontSize',  12)
        end
        plot(t,rateBaseHist,'LineStyle','--','Color',[.5,0.5,0.5])
        plot([0,0], [0, rmax], '--k'), hold off
        
        gf = figure(ifigRas);
        sUnits=num2str(un(1));
        for iu=2:numel(un)
            sUnits=sprintf('%s,%2d',sUnits,un(iu));
        end
        suptitle(sprintf('Mouse %s - session %3d - rec %s units %s',num2str(vp.mouse),str2num(vp.sess),vp.rec,sUnits));
        
        
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
        
%             save(fullfile(fn.fold_an_sess,sprintf('%s_warped_resp.mat',fileBase)), '-struct', 'vd','resp')
%             print(gf,'-depsc',fullfile(fullfile(fn.fold_an_sess,sprintf('%s_warped_resp.eps',fileBase))))
% 
        resp=vd.resp;
    end % ks 1:numel(stim)
end  %resp = view_rasters(sType,un)







%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function norm_respirationTimes = get_sniff_data(trial)

% Go through each trial and collect all respiration cycles, then find the
% mean sniff lengths for number of cycles that fit in window
allPhaseLengths=[];
for it=1:numel(trial)
    tr=trial(it);
    
    zeros = tr.sniffZeroTimes;
%     zeros = round(tr.sniffParabZeroTimes(1:2,:));
    if isempty(zeros) 
        continue
    end

    % Find ii, index of first inhalation after stim on
    if tr.odorConc >0
        % ii is index of first inhalation after final valve open
        ii = find(zeros(1,:)>round(tr.odorTimes(1)*1.009),1);
    elseif diff(tr.laserTimes)>0 && tr.odorConc==0
        ii = find(zeros(1,:)>=round(tr.laserTimes(1)*1.009),1);
    end  
    if isempty(ii)
        continue
    end

    % Get 2 sniffs before and 2 sniffs after stimulus
    respirationTimes = zeros(:,ii-2:ii+2);
    respirationTimes(2,end) = nan;
    
    % Check for possible missed cycles
    cyclePeriods = diff(respirationTimes(1,:));
    pauseLengths = respirationTimes(1,2:end) - respirationTimes(2,1:end-1);
    if any(cyclePeriods>1000) || any(pauseLengths>600)
        warning('sniff_analysis may have missed a cycle')
        continue
    end
    
    respirationTimes = reshape(respirationTimes,1,2*size(respirationTimes,2));
    
    phaseLengths = diff(respirationTimes);
    allPhaseLengths = [allPhaseLengths reshape(phaseLengths(1:end-1),2,(numel(phaseLengths)-1)/2)];
        
end

if size(allPhaseLengths,2)<numel(trial)
    error('Could not collect enough sniff cycles to get a good average')
end

mean_phaseLengths = round(mean(allPhaseLengths,2))

norm_respirationTimes = [-sum(mean_phaseLengths) -mean_phaseLengths(2)...
    0 mean_phaseLengths(1) sum(mean_phaseLengths)...
    sum(mean_phaseLengths)+mean_phaseLengths(1) 2*sum(mean_phaseLengths)]


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%




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