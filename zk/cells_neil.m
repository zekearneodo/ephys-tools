% This script goes through the cells database (/units) and picks up the
% metadata of units following one criterion (mouse, sessions, features of
% cells, through the conditions inn keepCells.
% With all those cells, it runs:
% - neil_trial_structure: picks up all the trials
% - merge_clusters: gets all the spikes of the clusters listed in the cluList
%   of the cell
% It saves:
% - 
function [cn, cellsArray]= cells_neil(doit)
global cn

    cn.units_meta    = @units_meta;
    cn.just_a_raster = @just_a_raster;
    cn.make_rasters  = @make_rasters;
    
    if nargin>0 && doit==1
        ffn=file_names();
        % load all cells in a cell structure
        mice = {'ZKawakeM72','KPawakeM72'};
        %mice = {'KPawakeM72'};
        %load all the cells into an array of cells
        disp(wrap_message('Getting cells and trial structures for Neil','*'));
        cellsArray = [];
        for is=1:numel(mice)
            cs=mice{is};
            cellFiles=dir(fullfile(ffn.fold_unit_db, [cs '*.mat']));
            allCells = arrayfun(@(x) load(fullfile(ffn.fold_unit_db,x.name)),cellFiles,'UniformOutput',false);
            cellsArray = [cellsArray [allCells{:}]];
        end
        
        %filter the cells by session (knowing when the experiments begun)
        keepCells =  (strcmpi('ZKawakeM72',{cellsArray.mouse}) & [cellsArray.sess]>3) | (strcmpi('KPawakeM72',{cellsArray.mouse}) & [cellsArray.sess]>5);
        
        cellsArray(~keepCells)=[];
        cellsArray(~([cellsArray.quality]==1))=[];
        % now you got an array of cells for the suffixes
        %with the cells selected, make all the units and place them ein the
        %export_data folder
        units_meta(cellsArray)
        % now go through all those cells and:
        % - find them in all the recs they appear in
        % - make the raster for every rec
        % - append it to the cell's structure
    end

end

% go to the trial structure and make a copy with the adecquate spikeTimes
function units_meta(cellsArray)
% - get the trial structures, empty the spikes and write it in the export
%   folder
% - get the sniffs and also place them in the export folder
% - get the clusters that make up a unit, merge them and create a unit file
% Input:
%   cellsArray : array of unit metadata structures (the output of
%                getUnits.py)
% Output:
%

mice = unique({cellsArray.mouse});
for im=1:numel(mice)
    mouse=mice{im};
    disp(wrap_message(['Mouse ' mouse],'.'));
    %get the sessions
    cellsMouse = cellsArray(strcmpi(mouse,{cellsArray.mouse}));
    sessList = sort(unique([cellsMouse.sess]'));
    for is=1:numel(sessList)
        sess=sessList(is);
        fprintf('Session %d:\n',sess)
        sessCells = cellsMouse([cellsMouse.sess]==sess);
        %select the recs in this sess that have cells
        recList = unique({sessCells.rec});
        for ir=1:numel(recList)
            rec=recList{ir};
            fprintf('   rec %s ... ',rec)
            %copy the sniff
            sn = load(fn.rsm_data,'Sniff');
            save(fullfile(fn.fold_exp_data,sprintf('%s_%03d_%s_sniff.mat',mouse,sess,rec)),'-struct','sn','Sniff');
            fn=file_names(mouse,sess,rec);
            load(fn.trial)
            trial=rmfield(trial,'spikeTimes');
            save(fn.exp_trial,'trial');
            %get all the units in the rec
            recCells = sessCells(strcmpi(rec,{sessCells.rec}));
            
            % go through all the units and make the unit (merge when multi-cluster)
            unit = [];
            for ic=1:numel(recCells)
                theCell  = recCells(ic);
                fn=file_names(theCell.mouse,theCell.sess,theCell.rec);
                oneUnit = merge_clusters(theCell.mouse,theCell.sess,theCell.rec,theCell.clu);
                oneUnit.uId = theCell.uId;
                oneUnit.sessCell = theCell.sessCell;
                unit = [unit oneUnit];
                save(fn.exp_spikes,'unit');
            end
            fprintf('(%d units)\n',numel(recCells));
        end
    end
    fprintf('\n');
end
save(fullfile(fn.fold_exp_data,'unitsmeta.mat'),'cellsArray');
end

function make_rasters(cellsArray)
% goes through an array of cells meta and get the rasters for all the recs
% it appears in
% Input:
%   cellsArray : array of unit metadata structures (the output of
%                getUnits.py)
%
cells_uId = unique({cellsArray.uId});

for ic = 1:numel(cells_uId)
    this_cell_instances = cellsArray(strcmpi(cells_uId(ic),{cellsArray.uId}));
    %for all the instances of this cell (recs it is in)
    %gather the rasters.
    raster = arrayfun(@(x) just_a_raster(x.mouse,x.sess,x.rec,x.sessCell),this_cell_instances);
    
    %data for the unit
    %quick check for screw ups in following the cell through recs
    %if the cell is litral its litral in all recs
    one_cell.light = unique([this_cell_instances.light]);
    if length(one_cell.light)>1
        warning('Mismatch in light responsiveness of cell %s trhough recs. Skipping cell.',cells_uId(ic))
    end
    
    one_cell.mouse  = this_cell_instances(1).mouse;
    one_cell.sess   = this_cell_instances(1).sess;
    one_cell.uid    = sprintf('%s_%03d_%03d',one_cell.mouse,one_cell.sess,this_cell_instances(1).sessCell);
    
    one_cell.raster = raster;

    %save it
    fn=file_names(one_cell.mouse,one_cell.sess);
    cellFn=fullfile(fn.fold_exp_data,sprintf('%s_cell.mat',one_cell.uid));
    save(cellFn,'-struct','one_cell');
end

end

function [raster] = just_a_raster(mouse,sess,rec,unitSessNumber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script for creating big raster of all the trials for a cell for neil,
%%% using the trial structure from export_data
%returns one raster for one particular rec in which the cell appears.
%unitId is the unit identifier trhough the session.

    %unit is a unit of the type of neil units
    %get the data of the unit, get all the trials the unit is in (from neil
    %trials, and make a raster centered on the first inhale after onset of
    %odor
    % 
    %unitNumber is which unit of that rec you want to get the raster.
    
    fn = file_names(mouse,sess,rec);
    trial = neil_trial_structure(mouse,sess,rec);
    save(fn.exp_trial);
    load(fn.exp_spikes);
    
    %meta data
    load(fullfile(fn.fold_exp_data,'unitsmeta.mat'));
    %look up the cell
    cellsOfRec =cellsArray( find( strcmpi(mouse,{cellsArray.mouse}) & strcmpi(rec,{cellsArray.rec}) & [cellsArray.sess]==sess));
    thisCell = cellsOfRec([cellsOfRec.sessCell]==unitSessNumber);
    unitNumber = find(strcmpi(thisCell.uId,{cellsOfRec.uId})); %number of cell amongst the cells of the rec
    thisUnit = unit(unitNumber);
    %     %select a particular odor, for debugging purposes
%     trial=trial(strcmpi('2-hydroxyacetophenone',{trial.odorName}));
%     nt = numel(trial);

    
    %the trials in exp are already refined (output of neil_trial_structure)
    %get that trial structure and go trial by trial adding rows to the big
    %raster
    t1 = -200;
    t2 = 2500;
    
    nt=numel(trial);
    odors   = {trial.odorName};
    concs   = [trial.odorConc];
    trialId = {trial.id};
    tVec = (t1:t2);
    spikes = zeros(nt,length(tVec));
    t0Vec = nan(1,numel(trial));
    sp = round(thisUnit.times);
    
    x    = zeros(1e5,1);    y    = zeros(1e5,1);
    nsp = 0;

    for it = 1:nt
        %get all the spikes for this unit within the window tVec
        ii = find(trial(it).sniffZeroTimes(1,:)>trial(it).odorTimes(1)*1.009,1);
        if isempty(ii) || ii<1
            continue
        end
        t0 = trial(it).sniffZeroTimes(1,ii) + trial(it).start;
        spikeTimes = sp((sp>t0+t1)&(sp<t0+t2))-t0;
        t0Vec(it) = t0;
        
        % compact raster
        n=numel(spikeTimes);
        x(nsp+(1:n)) = spikeTimes(:);
        y(nsp+(1:n)) = it*ones(1,n);
        nsp = nsp + n;
        if n>0
        spikes(it,spikeTimes-t1) = 1;   
        end
    end
    x = x(1:nsp);
    y = y(1:nsp);
    
    figure
    plot(x,y,'.','MarkerSize',7);
    
    %%%%%%%%%%%%%%%
    % Add cellId?  mouse_sess_rec_int --> int=unique integer of cell 
    %                                         in session from unitDb
    raster.odors   = odors;
    raster.concs   = concs;
    raster.trialId = trialId;
    raster.spikes  = spikes;
    raster.t       = tVec;
    raster.t0      = t0Vec;
    raster.rec     = rec;
    raster.cell    = thisCell;
    
    
    %for quick debugging of rasters
    raster.x = x;
    raster.y = y;
    
end

function trials = neil_trial_structure(mouse,sess,rec)
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


    