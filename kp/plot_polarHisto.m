function plot_polarHisto(mouse,sess,rec,cellId)
% Function will load in the data structures as formatted for Neil, 
fn=file_names(mouse,sess,rec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Load sniffs and warp spikes - BASELINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(fn.base_folder,'data_Neil',sprintf('%strialsBase.mat',fn.basename_an)));
load(fullfile(fn.base_folder,'data_Neil',sprintf('%sspikesBase_cell%i.mat',fn.basename_an,cellId)));

binSize=20;

% Collect all respiration cycles and find mean lengths
allPhaseLengths=[];
for itr=1:numel(spikesBase(1).trialId)
    trBase=trialsBase(itr);
    
    inhOnsets   = find(diff(trBase.sniffPhase)>1) +1;
    pauseOnsets = find(diff(trBase.sniffPhase)<-1) +1;
    if ~isempty(inhOnsets)
        fullCycles = find(pauseOnsets>inhOnsets(1),1):find(pauseOnsets<inhOnsets(end),1,'last');
        respirationTimes = [inhOnsets; pauseOnsets(fullCycles) nan];
        if min(find(inhOnsets(2:end)-pauseOnsets(fullCycles(1:end))>600))>1
            respirationTimes = reshape(respirationTimes,1,2*size(respirationTimes,2));
            trialsBase(itr).respTimes = respirationTimes;
            
            phaseLengths = diff(respirationTimes);
            allPhaseLengths = [allPhaseLengths reshape(phaseLengths(1:end-1),2,(numel(phaseLengths)-1)/2)];
        else
            trialsBase(itr).respTimes = [];
        end
    end
end

normPhaseLengths = round(mean(allPhaseLengths,2));
range_inhLen = [min(allPhaseLengths(1,:)) max(allPhaseLengths(1,:))];
range_pauseLen = [min(allPhaseLengths(2,:)) max(allPhaseLengths(2,:))];

figure(1) 
cdfplot(allPhaseLengths(1,:))

%% Find central chunk of inhalation times **note: might separate into 2+
% distributions!
%     inh_range = [70 100];
    inh_range = [85 105];
    sub_idxs = find(allPhaseLengths(1,:)>=inh_range(1) & allPhaseLengths(1,:)<=inh_range(2));
    normPhaseLengths = round(mean(allPhaseLengths(:,sub_idxs),2));


%%
% ** (debug by plotting rando trials === paste here from below)



% Now warp spike times 
iwp=0;
spikes_warped=zeros(1,sum(normPhaseLengths));
spikes_phaseH=[];
for itr=1:numel(spikesBase(1).trialId)
    respirationTimes = trialsBase(itr).respTimes;
    if numel(respirationTimes)<3
        continue
    end
    spikes_real = find(spikesBase.spikes(itr,:));
    
    for iph=1:2:numel(respirationTimes)-2
        phaseBounds = respirationTimes(iph:iph+2);
        if diff(phaseBounds(1:2))<inh_range(1) || diff(phaseBounds(1:2))>inh_range(2)
            continue
        end
        
        iwp=iwp+1; %%%% moved this outside of "isempty"
    	spikes_toWarp_inh   = spikes_real(spikes_real>phaseBounds(1) & spikes_real<phaseBounds(2));
        spikes_toWarp_pause = spikes_real(spikes_real>=phaseBounds(2) & spikes_real<phaseBounds(3));
        
        if ~isempty([spikes_toWarp_inh spikes_toWarp_pause])
        if ~isempty(spikes_toWarp_inh)
            warped_times = round((spikes_toWarp_inh-phaseBounds(1)).*(normPhaseLengths(1)/diff(phaseBounds(1:2))));
            spikes_warped(iwp,warped_times) = 1;  %%%%can be more than 1!!
        end
        if ~isempty(spikes_toWarp_pause)
            warped_times = round((spikes_toWarp_pause-phaseBounds(2)).*(normPhaseLengths(2)/diff(phaseBounds(2:3)))) + normPhaseLengths(1);
            spikes_warped(iwp,warped_times) = 1;
        end   
        spikes_phaseH = [spikes_phaseH spikes_warped(iwp,:)];
        end
    end 
end
rateHistBase = sum(spikes_warped,1)/iwp*1000;
extrams = rem(length(rateHistBase),binSize);
rateHistBase = sum(reshape(rateHistBase(1:end-extrams),binSize,length(rateHistBase(1:end-extrams))/binSize),1)/binSize;
% 
% figure(3); clf
% imagesc(spikes_warped)
% colormap('gray')
% 
% figure(4); clf
% subplot(2,1,1)
% polar(linspace(0,2*pi,length(rateHistBase)),rateHistBase)
% % polar(linspace(0,iwp*2*pi,iwp*sum(normPhaseLengths)),spikes_phaseH)
% subplot(2,1,2)
% plot(rateHistBase)
% ylim([min(rateHistBase)-2 max(rateHistBase)+2])
% % 
% 






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Load sniffs and warp spikes - ODORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cell=load(fullfile(fn.base_folder,'data_Neil',sprintf('%s_%s_%03i_cell.mat',fn.mouse,fn.sess,cellId)));
sniff=load(fullfile(fn.base_folder,'data_Neil',sprintf('%ssniff.mat',fn.basename_an)));
q=load(fullfile(fn.base_folder,'data_Neil',sprintf('%strial.mat',fn.basename_an)));
trial=q.trials; clear q;

spikes = cell.raster;
odors= unique(spikes.odors);
for io = 1:numel(odors)
    % Find trial indices of current odor
    odor_tr_idx = find(strcmpi(odors(io),spikes.odors));
    
    % Collect all respiration cycles and find mean lengths
    allPhaseLengths=[];
    for itr=odor_tr_idx
        tr=trial(itr);
        
        % Find t0 in order to identify sniffs during stimulus on
        stimON = find([tr.stim]==1,1);
        t0_inh1 = find(diff(tr.sniffPhase(stimON:end))>1,1)+stimON;
        
        if ~any(isnan(tr.sniffPhase))
        % Crop sniff data to match length of spike data collected
        sniffPhase = tr.sniffPhase(t0_inh1-1:end);

        inhOnsets   = find(diff(sniffPhase)>1) -1;
        pauseOnsets = find(diff(sniffPhase)<-1) -1;
        if numel(inhOnsets)>1 && ~isempty(pauseOnsets)
            firstFullCycle = [inhOnsets(1:2); pauseOnsets(find(pauseOnsets>inhOnsets(1),1)) nan];
            respirationTimes = reshape(firstFullCycle,1,2*size(firstFullCycle,2));
            trial(itr).respTimes = respirationTimes;
            
            phaseLengths = diff(respirationTimes);
            allPhaseLengths = [allPhaseLengths reshape(phaseLengths(1:end-1),2,(numel(phaseLengths)-1)/2)];      
        else
            warning('not enough sniff data to warp')
            trial(itr).respTimes = nan;
        end
        else
            itr
            warning('No sniffPhase data')
        end
    end
    
    normPhaseLengths = round(mean(allPhaseLengths,2));
    range_inhLen = [min(allPhaseLengths(1,:)) max(allPhaseLengths(1,:))];
    range_pauseLen = [min(allPhaseLengths(2,:)) max(allPhaseLengths(2,:))];
    
    figure(2)
    cdfplot(allPhaseLengths(1,:))
    
    %% Find central chunk of inhalation times **note: might separate into 2+
    % distributions!
    %     inh_range = [70 100];
    inh_range = [0 300];
    sub_idxs = find(allPhaseLengths(1,:)>=inh_range(1) & allPhaseLengths(1,:)<=inh_range(2));
    normPhaseLengths = round(mean(allPhaseLengths(:,sub_idxs),2));
    
    %%
    % Now warp spike times
    iwp=0;
    spikes_warped=zeros(1,sum(normPhaseLengths));
    for itr=odor_tr_idx
        spikes_real = find(spikes.spikes(itr,:))-200;
        
        respirationTimes = trial(itr).respTimes;
        if ~isempty(respirationTimes)
        for iph=1:2:numel(respirationTimes)-2
            phaseBounds = respirationTimes(iph:iph+2);
            if diff(phaseBounds(1:2))<inh_range(1) || diff(phaseBounds(1:2))>inh_range(2)
                continue
            end
            
            spikes_toWarp_inh   = spikes_real(spikes_real>=phaseBounds(1) & spikes_real<phaseBounds(2));
            spikes_toWarp_pause = spikes_real(spikes_real>=phaseBounds(2) & spikes_real<phaseBounds(3));
            
            if ~isempty([spikes_toWarp_inh spikes_toWarp_pause])
                iwp=iwp+1;
                if ~isempty(spikes_toWarp_inh)
                    warped_times = round((spikes_toWarp_inh-phaseBounds(1)).*(normPhaseLengths(1)/diff(phaseBounds(1:2)))) + 1;
                    spikes_warped(iwp,warped_times) = 1;
                end
                if ~isempty(spikes_toWarp_pause)
                    warped_times = round((spikes_toWarp_pause-phaseBounds(2)).*(normPhaseLengths(2)/diff(phaseBounds(2:3)))) + normPhaseLengths(1)+1;
                    spikes_warped(iwp,warped_times) = 1;
                end
            end
        end
        end
    end
    rateHist = sum(spikes_warped,1)/iwp*1000;

%     figure(3)
%     imagesc(spikes_warped)
%     colormap('gray')
    
    extrams = rem(length(rateHist),binSize);
    rateHist = sum(reshape(rateHist(1:end-extrams),binSize,length(rateHist(1:end-extrams))/binSize),1)/binSize;
    
    figure(3+io);
    suptitle(odors(io))
    subplot(2,1,1)
%     hold on
    polar(linspace(0,2*pi,length(rateHist)),rateHist,'r')
    hold on
    polar(linspace(0,2*pi,length(rateHistBase)),rateHistBase,'k')

    subplot(2,1,2)
    hold on
    plot(rateHistBase,'k')
    plot(rateHist,'r')
    ylim([min(rateHistBase)-2 max(rateHist)+2])
    
end






%% debug by plotting rando trials
% % trs = round(numel(spikesBase(1).trialId)*rand(1,50))
% % x=[];   y=[];    iy = 0;
% % for itr = trs
% %     if ~strcmp(spikesBase.trialId{itr},{trialsBase(itr).trialId})
% %         error('struct fields dont match')
% %     end
% %     spikes_real = find(spikesBase.spikes(itr,:))
% %     respirationTimes = trialsBase(itr).respTimes;
% %     if numel(respirationTimes)<3
% %         continue
% %     end
% %     iy=iy+1;
% %     respirationTimes = respirationTimes(1:3);
% %     
% %     spikes_inrange = spikes_real(spikes_real>(respirationTimes(1)-100) & spikes_real<(respirationTimes(3)+100))-respirationTimes(1)
% %     
% %     x = [x spikes_inrange];
% %     y = [y ones(1,numel(spikes_inrange)).*iy]
% % 
% % end
% % 
% % figure(2);
% % plot(x,y,'.k')
% % xlim([-200 800])

end






