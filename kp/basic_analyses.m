%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script for visualizing spikes aligned with the onset of stimuli

% outputs vr and response struct
% global vp contains the parameters for all the functions to call

%has functions:
% -visual_par_init
% -make_stimuli_set
% -view_rasters (default)
% -view_lfp (not called by default)

function [ba, resp] = basic_analyses(mouse,sess,rec,unit,sType)
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

ba.visual_par_init     = @visual_par_init;
ba.make_stimuli_set    = @make_stimuli_set;
ba.view_rasters        = @raster_psth;
ba.view_lfp            = @view_lfp;
ba.psths_oneOdor       = @psths_oneOdor;

if nargin < 5
    sType='laser'
end

vp = visual_par_init(mouse,sess,rec,sType);
vd = make_stimuli_set(vp);

if ~vp.surface && ~isempty(unit)
    resp = raster_psth(sType,unit);  %if spike data collected, do the raster visualization
end
end

function vp=visual_par_init(mouse,sess,rec,sType)
% visual_par_init(mouse,sess,rec,[unit])
%initializes or resets the parameters for the visualization
%that will be stored in the global vp;
%makes the stimuli set
global vd;

vp.warp = 0;  % 1 if want to plot in warped space, 0 if not
vp.sort = 0;  % 1 if want to sort trials by length of 1st inhalation, 0 if order presented


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


%========================
vp.bin = 20;    %========
%========================


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
    vp.par = {'odorName', 'odorConc', 'laserDur', 'laserPower', 'pulseOnsetDelay'};
elseif isfield(vp.tr,'pulseOffset')
    vp.tr_type = 2;                        % just KPawakeM72_004
    vp.par = {'odorName', 'odorConc', 'laserDur', 'laserPower', 'pulseGroup','pulseOffset'};
else
    vp.tr_type = 1;                        % old sessions
    vp.par = {'odorName', 'odorConc', 'laserDur', 'laserAmp'};
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


if vp.tr_type==2                 % for KPawakeM72 sess_004 only
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
stim = stim(2:end);

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


function resp = raster_psth(sType,un)
global vd vp;

% Get the stimuli presented
stim=vd.stim(vd.(sType).list);
sSort=vd.(sType).sort;

if isnumeric(stim(1).(sSort))
    [~, in_sort] = sort(-[stim.(sSort)]);
else
    [~, in_sort] = sort({stim.(sSort)});
end
stim = stim(in_sort);
ns=numel(stim);



%% =======================================================================
%   Plot raster and psth for requested units
%  =======================================================================

    tr  = vp.tr;
    bin = vp.bin;
    fn  = vp.fn;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%       Get mean respiration lengths        %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    snfsPreStim  = 1;
    snfsPostStim = 1;
    
    Sniffs = -snfsPreStim:snfsPostStim-1;
    
    [norm_respirationTimes,edges]  =  get_sniff_data(tr,Sniffs);
    normPhaseLengths = diff(norm_respirationTimes(1:3))

    t1  =  -snfsPreStim*sum(normPhaseLengths);
    t2  =  snfsPostStim*sum(normPhaseLengths);
    nt  =  t2-t1;
    tms =  t1:t2-1;  %time from t1:t2 in ms, not yet reshaped nor cropped

    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%            Get baseline rate              %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Averaging spikes in same size window (t1:t2) around inhalation onsets
    % for a few blank sniffs before each trial start.
%     bad_trials=[];
%     rateBase=zeros(1,nt);   spBase_byTrial=zeros(1,nt);
    nb=0; inhLenBase = []; timeToFirstInhOdor = []; laserTimeAfterInhDetected = [];
    nSnfs=3;
    
    nSniffs = numel(Sniffs);
    spikesBase_warped = zeros(1e4,nSniffs*sum(normPhaseLengths));
    x=[];     y=[];
    

    sniffLengths   = nan(numel(tr),nSnfs);
    inhLengthsStim = nan(numel(tr),nSnfs+1);
    for kt = 1:numel(tr)
        trial=tr(kt);
        if isempty(trial.sniffParabZeroTimes) || isempty(trial.spikeTimes) || isempty(trial.sniffZeroTimes)
            continue
        end
        sp = sort(round(vertcat(trial.spikeTimes{un})))';
%         if numel(sp>-500 & sp<1500)<3
%             aaa=573;
%             bad_trials = [bad_trials kt];
%             continue
%         end
                if trial.odorConc>0
                    iFirstSniff = find(trial.sniffZeroTimes(1,:)>trial.odorTimes(1)*1.009,1,'first');
                    t0Base = trial.sniffZeroTimes(1,iFirstSniff);  %time of first sniff after stimulus
                    
                    timeToFirstInhOdor = [timeToFirstInhOdor t0Base-trial.odorTimes(1)];
                    
                elseif trial.laserAmp>0 && any(strcmpi(trial.odorName,{'none','empty','dummy','not_an_event'}))
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

                % get durations of the three sniff cycles and inhalations
                %   after stimulus onset
                for nS=nSnfs-(0:nSnfs-1)
                    try
                        i_zeros = trial.sniffZeroTimes(1,iFirstSniff:iFirstSniff+nS);
                        p_zeros = round(trial.sniffParabZeroTimes(2,iFirstSniff:iFirstSniff+nS));
                        inhLengthsStim(kt,:) = p_zeros-i_zeros;
                        sniffLengths(kt,:)   = diff(i_zeros);
%                         Szeros = [i_zeros; p_zeros];
%                         sniffLengths(kt,1:nS) = diff(trial.sniffParabZeroTimes(1,iFirstSniff:iFirstSniff+nS));
%                         inhLengthsStim(kt,1:nS+1) = diff(trial.sniffParabZeroTimes(1:2,iFirstSniff:iFirstSniff+nS));
                        break
                    catch
                        continue
                    end
                end

                     
            % Now add up the spikes in as many non-overlapping windows as fit in the pre-stimulus period

                % Instead of finding windows by time (nt) get 2 or 3
                % cycles at a time.
                windows = floor((iFirstSniff-1)/nSniffs);
                
                for ib = 1:windows
%                    respirationTimes = trial.sniffParabZeroTimes(1:2,iFirstSniff-(nSniffs*ib):iFirstSniff-(nSniffs*(ib-1)))
                   respirationTimes = round([trial.sniffZeroTimes(1,iFirstSniff-(nSniffs*ib):iFirstSniff-(nSniffs*(ib-1))) ; trial.sniffParabZeroTimes(2,iFirstSniff-(nSniffs*ib):iFirstSniff-(nSniffs*(ib-1)))]);
                   
                   if any(diff(respirationTimes(1,:))>900)
                       continue
                   end
                   nb=nb+1; %keep track of the number of *segments* added to baseline

%                    ur = poissrnd(0.05.*ones(1,10000));
%                    ur = rand(1,10000);
%                    ur = poissrnd(sin([-4000:6000]*0.012)+4);
%                    clear sp
%                    sp=zeros(1,length(ur));
%                    sp(ur<=0.06) = 1;
%                    sp=find(sp>0)-4000;
                   
                   for iw=1:nSniffs
                       
                       phaseBounds = [respirationTimes(1,iw) respirationTimes(2,iw) respirationTimes(1,iw+1)];
                       if (phaseBounds(3)-phaseBounds(1))<edges(1) || (phaseBounds(3)-phaseBounds(1))>edges(2)
                           continue
                       end
                       
%                        clear sp_toWarp_inh sp_toWarp_pause
%                        spikes_toWarp_pause = [];
                       
                       sp_toWarp_inh   = sp(sp>phaseBounds(1) & sp<phaseBounds(2));

                           if ~isempty(sp_toWarp_inh)
                               clear warped_times
                               warped_times = round((sp_toWarp_inh-phaseBounds(1)).*(normPhaseLengths(1)/diff(phaseBounds(1:2))));
                               warped_times = warped_times + norm_respirationTimes(2*iw-1);
                               for iuwt=unique(warped_times)
                                   spikesBase_warped(nb,tms==iuwt) = sum(warped_times==iuwt);
                               end
                               x = [x warped_times];
                               y = [y nb.*ones(1,numel(warped_times))];
                           end
                           
                       sp_toWarp_pause = sp(sp>=phaseBounds(2) & sp<phaseBounds(3));
                       
                           if ~isempty(sp_toWarp_pause)
                               clear warped_times
                               warped_times = round((sp_toWarp_pause-phaseBounds(2)).*(normPhaseLengths(2)/diff(phaseBounds(2:3)))) + normPhaseLengths(1);
                               warped_times = warped_times + norm_respirationTimes(2*iw-1);
                               for iuwt=unique(warped_times)
                                   spikesBase_warped(nb,tms==iuwt) = sum(warped_times==iuwt);
                               end
                               x = [x warped_times];
                               y = [y nb.*ones(1,numel(warped_times))];
                           end
                         
                   end 
                end

%                 %%% Baseline for cumulative count latency measure
% %                 spBase_cC = spBase;
%                 %  spBase_cC = spBase(spBase>0);
% %                 for cc = 1:numel(spBase_cC)
% %                     cumulativeCountBase(nb,(spBase_cC(cc)-t1):nt) = cc;
% %                 end
    end
    spikesBase_warped = spikesBase_warped(1:nb,:);
%     figure; imagesc(-spikesBase_warped)
%     colormap('gray')
    
    rateBase_warped=sum(spikesBase_warped,1)/nb;

%     rateBase=rateBase/nb; %divide spk times collected by number of segments used

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
    

    %% =======================================================================
    %%  For each stimulus find spikes in analysis window for each trial
    %   Store data for raster (x,y) and histogram (rate, reshaped)
    
%     % Set up the figure    %1 fig, odors as panels
%     ifigRas=figure(1); clf
%     scrsz = get(0,'ScreenSize');
%     set(ifigRas,'Position',[1 (scrsz(4)/2) 3*scrsz(3)/4 scrsz(4)/2],...
%         'Nextplot','add');
%     nSubPlots=double(2*numel(stim)); hS=zeros(1,nSubPlots);
%     figRows=2;  figCols=ceil(nSubPlots/figRows);
%     for iSubFig=1:nSubPlots
%         hS(iSubFig)=subplot(figRows,figCols,iSubFig);
%         set(hS(iSubFig),'Nextplot','add');
%     end
    
    % separate fig each odor
    ifigRas=figure(4);
%     scrsz = get(0,'ScreenSize');
%     set(ifigRas,'Position',[1 (scrsz(4)/2) 3*scrsz(3)/4 scrsz(4)/2],...
%         'Nextplot','add');
    nSubPlots=2; hS=zeros(1,nSubPlots);
    figRows=2;  figCols=ceil(nSubPlots/figRows);
    for iSubFig=1:nSubPlots
        hS(iSubFig)=subplot(figRows,figCols,iSubFig);
        set(hS(iSubFig),'Nextplot','add');
    end
    
    odorcolors = {'b' 'k' [0.69 0.23 0.93] 'c' 'y' 'g' 'r'};
    
    ip  = 0;
    ktall=0;
        
    
    for ks = 2
        
        if stim(ks).laserDur>0 && stim(ks).odorConc>0
            continue
        end
        
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

        
        
        
        
        
        
%         rmax = round(3*max(rateBase_warped)*1000);
        rmax = 50;

        
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%     Now get stimulus-based properties     %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        rate = zeros(1,nt);

        %         cumulativeCount = zeros(1e3,nt);
        x = [];   y = [];  spikes_warped = zeros(1,nSniffs*sum(normPhaseLengths));
        kt   = 0;
        
        ISI_raw = [];
        spkCount=[];	
        firstSniff_stim=[];
        
        for it = [stim(ks).in_tr]
            if isempty(tr(it).start) || isempty(tr(it).sniffParabZeroTimes) || isempty(tr(it).sniffZeroTimes)
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
                
            elseif strcmp(sType,'laser')
                t0 = round(tr(it).laserTimes(1)*1.009);
                kt = kt + 1;
            end
            

            
            % Instead of finding windows by time (nt) get 2 or 3
            % cycles at a time.
            try
%                 respirationTimes = round(tr(it).sniffParabZeroTimes(1:2,ii+(-snfsPreStim:snfsPostStim)));
                respirationTimes = round([tr(it).sniffZeroTimes(1,ii+(-snfsPreStim:snfsPostStim)) ; tr(it).sniffParabZeroTimes(2,ii+(-snfsPreStim:snfsPostStim))]);
            catch
                continue
            end

            if any(diff(respirationTimes(1,:))<edges(1)) || any(diff(respirationTimes(1,:))>edges(2))
                continue
            end
                
            kt = kt + 1;
%             timetoFirstInh(kt) = t0-tr(it).odorTimes(1)*1.009;
            firstSniff_stim(kt,:) = [diff(respirationTimes(:,snfsPreStim+1)) respirationTimes(1,snfsPreStim+2)-respirationTimes(2,snfsPreStim+1)];
            
% %             % Simulate spike data to check warping
% %             ur = poissrnd(0.05.*ones(1,10000));
% % %             ur = poissrnd(sin([-4000:6000]*0.012)+4);
%             ur = 5*kt:50:10000;
%             sp=zeros(1,1e4);
%             sp(ur) = 1;
%             sp=find(sp)-2000;
            
            for iw=1:nSniffs

                phaseBounds = [respirationTimes(1,iw) respirationTimes(2,iw) respirationTimes(1,iw+1)];
                if (phaseBounds(3)-phaseBounds(1))<edges(1) || (phaseBounds(3)-phaseBounds(1))>edges(2)
                    break
                end
                
                sp_toWarp_inh   = sp(sp>phaseBounds(1) & sp<phaseBounds(2));
                sp_toWarp_pause = sp(sp>=phaseBounds(2) & sp<phaseBounds(3));
                if ~isempty([sp_toWarp_inh sp_toWarp_pause])
                if ~isempty(sp_toWarp_inh)
                    clear warped_times
                    warped_times = round((sp_toWarp_inh-phaseBounds(1)).*(normPhaseLengths(1)/diff(phaseBounds(1:2))));
                    warped_times = warped_times + norm_respirationTimes(2*iw-1);
                    for iuwt=unique(warped_times)
                        spikes_warped(kt,tms==iuwt) = sum(warped_times==iuwt);
                    end
                    x = [x warped_times];
                    y = [y kt*ones(1,numel(warped_times))];
                end
                if ~isempty(sp_toWarp_pause)
                    clear warped_times
                    warped_times = round((sp_toWarp_pause-phaseBounds(2)).*(normPhaseLengths(2)/diff(phaseBounds(2:3)))) + normPhaseLengths(1);
                    warped_times = warped_times + norm_respirationTimes(2*iw-1);
                    for iuwt=unique(warped_times)
                        spikes_warped(kt,tms==iuwt) = sum(warped_times==iuwt);
                    end
                    x = [x warped_times];
                    y = [y kt*ones(1,numel(warped_times))];
                end
                tr_sp = find(spikes_warped(kt,:)>0)+t1;
                end
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%   Calculate ISI for each sp in trial
            ISI_raw(kt,:) = nan(1,nt);  %%%use a sparse matrix for this?
            if numel(tr_sp) > 1
                ISI_raw(kt,tr_sp(2:end)-t1) = diff(tr_sp);  % vector of ISIs
                ISI_raw(kt,tr_sp(1)-t1) = inf;
                quickind = find(~isnan(ISI_raw(kt,:)));
                quicklook = ISI_raw(kt,quickind);
            else
                quickind = [];
                quicklook = [];
            end
        end

        spkCount = mean(spkCount);
        
        rate = sum(spikes_warped,1)/kt;
%         rate_hz = rate*1000;

%         if kt<15
%             continue
%         end
        

        %         figure;
        %         cdfplot(timetoFirstInh)
        
        % Plot distributions of lengths of:
        if strcmp(sType,'odor')
        figure(3); hold on
        subplot(2,3,4);
        hist(sum(firstSniff_stim,2),100,'r') % whole cycle
        xlim([0 800])
        subplot(2,3,5); 
        hist(firstSniff_stim(:,1),100,'r') % inhalation period
        xlim([0 300])
        subplot(2,3,6); 
        hist(firstSniff_stim(:,2),100,'r') % pause+exhalation period
        xlim([0 800])
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%       Calculate response properties       %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Crop extra 1ms bins so time and rate vectors can be reshaped
        if ip==0
            extrams = rem(nt,bin);
            nt  =  t2-t1-extrams;
            tms =  t1:t2-1-extrams;  %time from t1:t2 in ms, not reshaped
            t   =  mean(reshape(tms,bin ,nt/bin),1);
            
            vp.smooth.wsize = round(nt/8);   %window size for gaussian smoothing of histo for plotting
            vp.smooth.cutoff = 100;   %cutoff for the gaussian smoothing
            vp.smooth.stdev = 19358/(2*pi*vp.smooth.cutoff); %std for gaussian smoothing of histo for plotting
            smooth = vp.smooth;
        end
    
        
        %%%%    BASELINE    %%%%
        rateBaseHist = sum(reshape(rateBase_warped(1:end-extrams), bin, nt/bin),1)/(bin/1000); %units is Hz
        rateBaseHist_smooth = smoothts(rateBaseHist,'g',smooth.wsize,smooth.stdev);
        avgRateBase  = mean(rateBase_warped*1000) %units is Hz
        stdRateBase  = std(rateBase_warped*1000);
        rateBase_std_bin = std(reshape(rateBase_warped(1:end-extrams)*1000, bin, nt/bin),1);
        
        spkCountBase   = sum(rateBase_warped(tms>0 & tms<responseWindow));
        nSpikesBaseAvg = avgRateBase*responseWindow/1000;
        
        
        
        %%%%    STIMULUS    %%%%
        rateHist = sum(reshape(rate(1:end-extrams), bin, nt/bin),1)/(bin/1000); %units is Hz
        rateHist_smooth = smoothts(rateHist,'g',smooth.wsize,smooth.stdev);

        
        rate_subtrBase = rateHist-rateBaseHist;
      
        
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
        %%%    Interpolate point where rateHist crosses time-varying threshold
        
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
        
        

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%  Find first peak and max peak of the binned histogram
        
        % Manual find peaks code, so that threshold can vary with time
        t_stim=t(t>0);                  rateBase_stim = rateBaseHist(t>0);
        rateHist_stim=rateHist(t>0);    rateBase_std_stim=rateBase_std_bin(t>0);
        
%         manPks_allin = 1+find(diff(sign(diff(rateHist_stim)))==-2);
%         manPks = rateHist_stim(manPks_in) > (1.5*rateBase_std_stim(manPks_in) +rateBase_stim(manPks_in));
% Same as Shusterman: peaks with rate 2x baseline at that tbin
        [~,pkInd]=findpeaks(rateHist_stim,'MINPEAKHEIGHT',min(rateBase_stim),'SORTSTR','none');
        manPks = rateHist_stim(pkInd) >= (2*rateBase_stim(pkInd));
        manPks_in = pkInd(manPks);
        
        manPks_time = t_stim(manPks_in);
        manPks_rate = rateHist_stim(manPks_in);
        
        
        if isempty(manPks_rate)
            firstPk_time=nan
            firstPk_rate=nan;
            maxPk_time=nan
            maxPk_rate=nan;
        else
            firstPk_time = round(manPks_time(1))
            firstPk_rate = manPks_rate(1);
            maxPk_time = round(t_stim(find(rateHist_stim==max(manPks_rate),1)))
            maxPk_rate = max(manPks_rate);
        end
        
        % Count number of spikes in first sniff with stim
        nSpikes=round(sum(rate(tms>=0 & tms<=sum(normPhaseLengths))));
        
        
        %-----------------------------------------------------------
        %%%%     SHUSTERMAN/RINBERG precise latency calc
        ISItoplot = []; nsp = 0; firstISI = nan(kt,1); nSpkISI = zeros(kt,1);

        % If there was a sharp event during the first sniff
        if ~isnan(firstPk_time) && firstPk_time<=sum(normPhaseLengths)
            
            % Define isi threshold for this cell-odor pair
            fm = firstPk_rate;
            isi_thr = 1.5/fm*1000; %ms
            
            % select time window around response peak
            latencyWin = [firstPk_time-30 firstPk_time+60];
            latencyWin(latencyWin<0) = 0;  %correct windows beginning before t0
            
            for jt = 1:size(spikes_warped,1)
                tr_sp = find(spikes_warped(jt,:))+t1;
                tr_sp = tr_sp(find(tr_sp>=latencyWin(1),1):min(length(tr_sp),find(tr_sp<=latencyWin(2),1,'last')+1));
                trISI = diff(tr_sp);
                if ~isempty(find(trISI<isi_thr,1))
                    %time of first spike where next follows by less than threshold ms
                    firstISI(jt) = tr_sp(find(trISI<isi_thr,1));
                end
            end
        end
        
        latencyISI = mean(firstISI(isfinite(firstISI)))
        jitter = std(firstISI(isfinite(firstISI)))
        reliability = numel(firstISI(isfinite(firstISI)))/kt
        
%         if strcmp(sType,'odor')
%             firstISI_sorted = firstISI(in_inh1);
%         end
        
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%   Save properties to response structure   %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        vd.resp(ks).rateBase_warped     = rateBase_warped;
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
        vd.resp(ks).reliabilityISI      = reliability;
        vd.resp(ks).latency             = latency_interpolated;
        vd.resp(ks).inhibitionDetected  = inhibition_interpolated;
        %
        vd.resp(ks).normPhaseLengths    = normPhaseLengths;
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
        gs = hS(1);
        subplot(gs); hold on
%         if strcmp(sType,'odor')
%         for isn = 1:2:2*nSniffs
%             inh_plot_x = [norm_respirationTimes(isn).*ones(1,2) norm_respirationTimes(isn+1).*ones(1,2)];
%             inh_plot_y = [0 kt+1 kt+1 0];
%             fill(inh_plot_x,inh_plot_y,[0.95 0.95 0.95],'LineStyle',':','LineWidth',0.001)
%         end
%         end
        box off
        plot(x,y+ktall, '.', 'Color',odorcolors{ip}, 'MarkerSize',6)
        set(gs, 'XLim', [t1 t2],...
            'YLim',      [0 kt+1+ktall],...
            'YTick',     [0 kt+ktall],...
            'XTick',     norm_respirationTimes(find(mod(1:numel(norm_respirationTimes),2))),...
            'FontSize',  10);
%         plot(firstISI,1:length(firstISI), 'r.', 'MarkerSize',8)
        plot([0,0], [0, kt+ktall+1], '--k'), hold off
%         title(gs,stim_str,'FontSize', 14)
        stim_handle{ip} = stim_str;
        
        ktall = ktall + kt;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot histogram in lower subplot
%         gs = hS(ip+numel(stim));
        gs = hS(2);
        subplot(gs); hold on
        if ip==1
            br=area(t,rateBaseHist,'FaceColor',[.8,0.8,0.8],'LineStyle','none')
        end
        rh(ip)=plot(t, rateHist, ...
            'LineWidth', 3, ...
            'Color',     odorcolors{ip}, ...
            'LineStyle', '-')
        hold on
        box off
        set(gs,'XLim',       [t1, t2], ...
            'XTick',         norm_respirationTimes, ...
            'YLim',          [0 rmax], ...
            'YTick',         0:50:rmax, ...
            'FontSize',      12);
%         plot(firstPk_time,firstPk_rate,'c*')
%         plot(latencyISI,0.75*firstPk_rate,'r*','MarkerSize',8)
% (find(mod(1:numel(norm_respirationTimes),2)))

%         set(br,'Color',[.4,0.4,0.4])
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
        
%         vd.resp(ks).visInspect   =  0;
%         save(fullfile(fn.fold_an_sess,sprintf('%s_warped_resp.mat',fileBase)), '-struct', 'vd','resp')
%         print(gf,'-depsc',fullfile(fullfile(fn.fold_an_sess,sprintf('%s_warped_resp.eps',fileBase))))
%         
        resp=vd.resp;
    end % ks 1:numel(stim)
    legend(rh,stim_handle)
% figure; hold on
% plot([resp.reliabilityISI],'k.--')
% plot(([resp.jitterISI]-min([resp.jitterISI]))./(max([resp.jitterISI]-min([resp.jitterISI]))),'b.--')
% ylim([0 1])
end  %resp = view_rasters(sType,un)







%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [norm_respirationTimes,edges] = get_sniff_data(trial,Sniffs)

% Go through each trial and collect all respiration cycles, then find the
% mean sniff lengths for number of cycles that fit in window
nC = 0;
for it=1:numel(trial)
    tr=trial(it);
    if isempty(tr.start)
        continue
    end
    if isempty(tr.sniffParabZeroTimes) 
        if it==numel(trial)
            warning('No parab zeros saved in trial struct for this mouse! Run sniff_analysis again.')
        end
        continue
    end
    
%     Szeros = round(tr.sniffParabZeroTimes(1:2,:));
    
    i_zeros = tr.sniffZeroTimes;
    p_zeros = round(tr.sniffParabZeroTimes(1:2,:));
    try
        Szeros = [i_zeros(1,:); p_zeros(2,:)];
    catch
        warning('zero crossings and parab fits dont have same number of elements')
    end

    % Get lengths from each sniff cycle in this trial
    for icycle = 1:size(Szeros,2)-1
        thisCycle   = Szeros(:,icycle:icycle+1);
        cyclePeriod = diff(thisCycle(1,:));
        if cyclePeriod>1000
            continue
        end
        nC = nC +1;
        phaseLengths(nC,1) = diff(thisCycle(:,1));
        phaseLengths(nC,2) = thisCycle(1,2) - thisCycle(2,1);
    end
end

% Keep only 10th to 90th percentile of whole cycle lengths
edges = prctile(sum(phaseLengths,2),[5 95]);
phaseLengths = phaseLengths(sum(phaseLengths,2)>edges(1)&sum(phaseLengths,2)<edges(2),:);

% Plot distributions of lengths of:
figure(3); clf
subplot(2,3,1)
hist(sum(phaseLengths,2),100) % whole cycle
xlim([0 1000])
title('whole sniff lengths')
subplot(2,3,2)
hist(phaseLengths(:,1),100) % inhalation period
xlim([0 300])
title('inhalation lengths')
subplot(2,3,3)
hist(phaseLengths(:,2),100) % pause+exhalation period
xlim([0 800])
title('pause/exh lengths')


% Calculate means for each phase to use for warping spike data
mean_phaseLengths = round(mean(phaseLengths,1));
norm_respirationTimes = [0 cumsum(reshape(mean_phaseLengths'*ones(1,(numel(Sniffs))),1,(numel(Sniffs))*2))] + Sniffs(1)*sum(mean_phaseLengths);

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
t1=-500;
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
            ii = find(tr(it).sniffParabZeroTimes(1,:)>tr(it).odorTimes(1),1);
            if isempty(ii) | ii<1
                continue
            end
            t0 = tr(it).sniffParabZeroTimes(1,ii);
            
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










function allResp = psths_oneOdor()
global rp
fn=file_names();

if ~isfield(rp,'qualityCells')
    %% Get all cells with quality>0
    nCells=0;
    find_qCells=[];
    cellBaseName='*.mat';
    cellsList=dir(fullfile(fn.fold_unit_db,cellBaseName));
    find_qCells=cellfun(@(x) getCell(fullfile(fn.fold_unit_db,x)),{cellsList.name},'UniformOutput',false);
    find_qCells(cellfun('isempty',find_qCells))=[];
    qCells(nCells+1:nCells+numel(find_qCells))=[find_qCells{:}];
    rp.qualityCells=qCells;
elseif isfield(vp,'qualityCells') && ~isempty(rp.qualityCells)
    qCells = rp.qualityCells;
end
nCells=numel(qCells);


q=load(fullfile(fn.base_folder,'stimuli','odorDb'));
odorDb=q.odor; clear q








end













    
function odorInfo=get_odor_info(voPar)
odorInfo=struct('name','','flows',[],'dillution',1,'vialConc',{},'vialId',{},'vial',{});

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
