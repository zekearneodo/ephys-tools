function [spikesBase,qCells] = assemble_baseline_spikeStruct(qCells)
%% This function creates a file of cell/spike data formatted for Neil. 
fn=file_names();

% Specify cell base name if you want to
%           need to generalize to optional specification


% Find cell meta data for the unit(s) you want to get spikes for. 
% If the cells you want are already a variable in the workspace, enter 
% qCells as an input argument.
if nargin<1
%     mouse = 'ZKawakeM72';
%     sess  = '004';
%     rec   = 'c';
    
    nCells=0;
    find_qCells=[];
    cellBaseName=[sprintf('%s_%s_%s_',mouse,sess,rec) '*.mat'];
%     cellBaseName=['*.mat'];
    cellsList=dir(fullfile(fn.fold_unit_db,cellBaseName));
    find_qCells=cellfun(@(x) getCell(fullfile(fn.fold_unit_db,x)),{cellsList.name},'UniformOutput',false);
    find_qCells(cellfun('isempty',find_qCells))=[];
    qCells(nCells+1:nCells+numel(find_qCells))=[find_qCells{:}];
end
nCells=numel(qCells);

% 
for ic = 1:nCells
    [~,remain] = strtok(qCells(ic).Id,'_');
    [sess,~] = strtok(remain,'_');
    fn=file_names(qCells(ic).mouse,sess,qCells(ic).rec);
    
    % Skip recs that had old method of aligning trials
    if any(strcmpi(fn.basename_an,{'KPawakeM72_006_a_' 'ZKawakeM72_001_c_'}))  
        continue
    end
    
    % Load trial struct for current cell
    q = load(fn.trial);
    trial=q.trial; clear q
    resp=qCells(ic).resp;
    
    units = [qCells(ic).clu]';   % this is cluster number(s) in rec
    
    % Try to load trialsBase struct (format for Neil) for the
    % mouse/sess/rec the current cell belongs to. 
    % If it doesn't exist, call function to create it.
    try
        q = load(fullfile(fn.base_folder,'data_Neil',sprintf('%strialsBase.mat',fn.basename_an)));
        trialsBase = q.trialsBase; clear q
    catch
        fprintf('\tRunning assemble_baseline_trialStruct for %s_%s_%s\n',qCells(ic).mouse,sess,qCells(ic).rec)
        [trialsBase] = assemble_baseline_trialStruct(qCells(ic).mouse,sess,qCells(ic).rec);
    end
    
    %% Extract spikes and save structure
    spikes = zeros(numel(trialsBase),numel(trialsBase(1).sniffPhase));
%     spikes_b = zeros(numel(trialsBase),numel(trialsBase(1).sniffPhase));
    otr=1;  tr_uids=[]; snPhase = [];
    for itr = 1:numel(trial)
        if otr<=numel(trialsBase)
            trialUId = [fn.basename_an 'trial' num2str(trial(itr).start)];
            if ~strcmp(trialUId,trialsBase(otr).trialId)
                continue
            else
                if ~isempty(trial(itr).spikeTimes)
                    allspikes = sort(round(vertcat(trial(itr).spikeTimes{units})));
%                         + trialsBase(otr).start;
%                     allspikes_b = sort(round(vertcat(trial(itr).spikeTimes{units})))...
%                         + trialsBase(otr).start;
%                     spikeTimes_b = allspikes_b(allspikes_b>trialsBase(otr).start...
%                         & allspikes_b<=(trialsBase(otr).start+numel(trialsBase(otr).sniffPhase)))...
%                         - trialsBase(otr).start;
                    t1fromTrPin = trialsBase(otr).start - trial(itr).start;
                    spikeTimes = allspikes(allspikes>t1fromTrPin...
                        & allspikes<=(t1fromTrPin+numel(trialsBase(otr).sniffPhase)))...
                        -t1fromTrPin;
                    spikes(otr,spikeTimes) = 1;
%                     spikes_b(otr,spikeTimes_b) = 1;
                    if ~any(isnan(trialsBase(otr).sniffPhase))
                        snPhase = [snPhase; trialsBase(otr).sniffPhase];
                    end
                    otr=otr+1;
                else
                    spikes(otr,:) = nan;
%                     spikes_b(otr,:) = nan;
                    otr=otr+1;
                end
            end
            tr_uids = [tr_uids; {trialUId}];
        else
            break
        end
    end
    if size(spikes,1) ~= numel(trialsBase)
        warning('didnt catch all odor trials?')
    else
        cellId = sprintf('%s_%s_%s_%i',qCells(ic).mouse,sess,qCells(ic).rec,qCells(ic).sessCell);
        allSpikesBase(ic).cellId = cellId;
        allSpikesBase(ic).spikes = spikes;
        allSpikesBase(ic).trialId = tr_uids;
        
%     figure; 
%     subplot(2,1,1);
%     imagesc(spikes); colormap('gray')
%     subplot(2,1,2);
%     imagesc(snPhase); colormap('gray')
    end
    
%     if strcmp(cellId,'KPawakeM72_016_a_16')
        spikesBase = allSpikesBase(ic);
%         save(fullfile(fn.base_folder,'data_Neil',sprintf('%sspikesBase_cell%i.mat',fn.basename_an,qCells(ic).sessCell)), 'spikesBase')
%     end
end
% Save all cells in one struct
% save(fullfile(fn.base_folder,'data_Neil',sprintf('spikesBase_lightrals_%s.mat',date)), 'allSpikesBase')

   


%_______________________________________________________________________%
    function theCell=getCell(unit_filename)
        theCell=load(unit_filename);
        if theCell.quality==1 && theCell.light==1 && theCell.odor==1
            theCell.resp = get_resp_struct(unit_filename,'odor');
        else
            theCell='';
        end
    end
end

