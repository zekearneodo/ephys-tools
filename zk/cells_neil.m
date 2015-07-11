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

%doit: 0 (or empty) [default] Don't run the rasters
%      1 make the rasters for neil (no noStimSniffs, for instance)
%      2 make the rasters for continuing analysis in python
global cn

cn.units_meta      = @units_meta;
cn.do_raster       = @do_raster;
cn.just_a_raster   = @just_a_raster_set;
cn.just_a_baseline = @just_a_baseline;
cn.make_rasters    = @make_rasters;
cn.neil_trial_structure = @neil_trial_structure;

cn.assemble_baseline_trial_tructure = @assemble_baseline_trial_tructure;
cn.assemble_baseline_trial_tructure_zk = @assemble_baseline_trial_tructure_zk;
cn.assemble_baseline = @assemble_baseline;

cn.no_stim_sniffs_raster = @no_stim_sniffs_raster;

if nargin==0
    doit = 0;
end

if doit>0
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
        cellsArray = [cellsArray [allCells{:}]]; %#ok<AGROW>
    end
    
    keepCells1 = zeros(1,numel(cellsArray));
    keepCells2 = zeros(1,numel(cellsArray));
    %filter the cells by session (knowing when the experiments begun)
    %affinity:
    keepCells1 =  (strcmpi('ZKawakeM72',{cellsArray.mouse}) & ([cellsArray.sess]>3 & [cellsArray.sess]< 15 | [cellsArray.sess]==914)...
        | strcmpi('KPawakeM72',{cellsArray.mouse}) & [cellsArray.sess]>13 & [cellsArray.sess]<17);

    %concentration:
    keepCells2 =  (strcmpi('ZKawakeM72',{cellsArray.mouse}) & [cellsArray.sess]>17 & [cellsArray.sess]< 28 ...
        | strcmpi('KPawakeM72',{cellsArray.mouse}) & [cellsArray.sess]>17);
    
    keepCells = keepCells1 + keepCells2;
    
    %for neil new emission (7/8/15)
    keepCells3 = strcmpi('KPawakeM72',{cellsArray.mouse}) & ([cellsArray.sess]==23 | [cellsArray.sess]==24);
    keepCells4 = strcmpi('ZKawakeM72',{cellsArray.mouse}) & ([cellsArray.sess]==28 | [cellsArray.sess]==29 | [cellsArray.sess]==15);
    
    keepCells = keepCells1 + keepCells2 + keepCells3 + keepCells4;
   
    
    if doit ==2
        %for neil debugging:
        keepCells = strcmpi('ZKawakeM72',{cellsArray.mouse}) & ([cellsArray.sess]==13);
    end
    
    cellsArray(~keepCells)=[];
    cellsArray(~([cellsArray.quality]==1))=[];
    %cellsArray(~([cellsArray.light]==1))=[];
    % now you got an array of cells for the suffixes
    %with the cells selected, make all the units and place them ein the
    %export_data folder
    units_meta(cellsArray);
    
    % now go through all those cells and:
    % - find them in all the recs they appear in
    % - make the raster for every rec
    % - append it to the cell's structure
    
    %cellsArray(~(strcmpi('ZKawakeM72',{cellsArray.mouse}) & [cellsArray.sess]<20))=[];
    make_rasters(cellsArray, doit);
    %save all the meta,
    save(fullfile(ffn.fold_exp_data,'cells_meta.mat'), 'cellsArray');
end

end

% go to the trial structure and make a copy with the adequate spikeTimes
function units_meta(cellsArray)
% - get the trial structures, empty the spikes and write it in the export
%   folder
% - get the sniffs and also place them in the export folder
% - get the clusters that make up a unit, merge them and create a unit file
% Input:
%   cellsArray : array of unit metadata structures (the output of
%                getUnits.py)
% Output:
% Writes files:
% - basename_sniff.mat      : sniff trace for the rec
% - basename_trial.mat      : trial structure for the rec
% - basename_trialsBase.mat : trial structure for the baseline
% - fn.exp_spikes           : spikes for the unit (one per unit)
% - unitsmeta.mat           : the array of units

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
            fn = file_names(mouse,sess,rec);
            sn = load(fn.rsm_data,'Sniff');
            save(fullfile(fn.fold_exp_data,sprintf('%s_%03d_%s_sniff.mat',mouse,sess,rec)),'-struct','sn','Sniff');
            fn=file_names(mouse,sess,rec);
            %do the trial structure for the baseline (using kristina's
            %program)
            %trialsBase = assemble_baseline_trial_tructure_zk(mouse,sess,rec);
            %save(fullfile(fn.fold_exp_data,sprintf('%strialsBase.mat',fn.basename_an)), 'trialsBase')
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

function make_rasters(cellsArray, doit)
% goes through an array of cells meta and get the rasters for all the recs
% it appears in
% Input:
%   cellsArray : array of unit metadata structures (the output of
%                getUnits.py)
%   doit       : 1 (or empty) make the rasters for neil
%                2 make them for cont. analysis
%cellsArray(~strcmpi('ZKawakeM72_010_001',{cellsArray.uId}))=[];

cp = cell_passport_tools;

if nargin==0
    doit = 1;
end

cells_uId = unique({cellsArray.uId});

for ic = 1:numel(cells_uId)
    this_cell_instances = cellsArray(strcmpi(cells_uId(ic),{cellsArray.uId}));
    %for all the instances of this cell (recs it is in)
    %gather the rasters.

    spikesBase = arrayfun(@(x) just_a_baseline(x.mouse,x.sess,x.rec,x.sessCell),this_cell_instances);
    %data for the unit
    %quick check for screw ups in following the cell through recs
    %if the cell is litral its litral in all recs
    clear one_cell;
    one_cell.light = unique([this_cell_instances.light]);
    if length(one_cell.light)>1
        warning('Mismatch in light responsiveness of cell %s trhough recs. Skipping cell.',cells_uId(ic))
    end
    
    one_cell.mouse  = this_cell_instances(1).mouse;
    one_cell.sess   = this_cell_instances(1).sess;
    one_cell.uid    = sprintf('%s_%03d_%03d',one_cell.mouse,one_cell.sess,this_cell_instances(1).sessCell);
    
    %get info of the recording to add to the cell
    sess_info = cp.load_info(this_cell_instances(1));
    rec = this_cell_instances(1).rec;
    rec_info = sess_info(strcmpi(rec, {sess_info.rec}));
    if isfield(rec_info, 'site_side')
        one_cell.is_ipsi = isempty(strcmpi(rec_info.site_side, 'contra'));
    else
        one_cell.is_ipsi = 0;
    end
    
    % add the rasters of each response
    %fields to store in the raster
    %field for storing -> name of field in read raster
    raster_fields.odor       = 'raster';
    raster_fields.light      = 'light_raster';
    raster_fields.odor_light = 'odor_light_raster';
    

    for i = 1:numel(this_cell_instances)
        x = this_cell_instances(i);
        raster = just_a_raster_set(x.mouse,x.sess,x.rec,x.sessCell, doit);
        if isempty(raster)
            continue
        end
        
        rf = fields(raster);
        for ifr = 1:numel(rf)
            read_field  = rf{ifr};
            write_field = raster_fields.(read_field);
            
            if ~isfield(one_cell, write_field)
                one_cell.(write_field) = raster.(read_field);
            else
                one_cell.(write_field) = [one_cell.(write_field) raster.(read_field)];
            end
        end
        
        %make the no stim sniff structure and save it
        if doit == 2
            noStimSniff = no_stim_sniffs_raster(x);
            fn=file_names(x.mouse,x.sess,x.rec);
            noStimFn = fullfile(fn.fold_exp_data,sprintf('%snoStimSniff.mat',fn.basename_an));
            save(noStimFn, 'noStimSniff');
        end

    end
    
    %save it
    fn=file_names(one_cell.mouse,one_cell.sess);
    cellFn=fullfile(fn.fold_exp_data,sprintf('%s_cell.mat',one_cell.uid));
    save(cellFn,'-struct','one_cell');
    
    %save the baseline
    baseFn = fullfile(fn.fold_exp_data,sprintf('%s_spikesBase.mat',one_cell.uid));
    save(baseFn,'spikesBase');
    
end

end

%create one raster from trial structure trial of events in sp_times, with
%stimulus type sType

function [raster] = do_raster(trial, sp, s_type)
%trial : trial structure
%sp : vector of spike times of the unit
%sType : stimulus type
t1 = -200;
t2 = 2500;

%before continuing, reduce the trials to only the ones that have the
%stimulus
align = [];
switch lower(s_type)
    % do only odor raster
    case 'odor'
        trial([trial.laserAmp]>0) = [];
        trial(strcmpi('none',{trial.odorName}) | strcmpi('empty',{trial.odorName}) | strcmpi('empty',{trial.odorName})| strcmpi('not_an_event',{trial.odorName})) = [];
        keep.odors      = {trial.odorName};
        keep.concs      = [trial.odorConc];
        align = 'sniff';
        %do only laser raster
    case {'laser', 'light'}
        trial([trial.laserAmp]<1) = [];
        trial( ~(strcmpi('none',{trial.odorName}) | strcmpi('empty',{trial.odorName})...
            | strcmpi('empty',{trial.odorName}) | strcmpi('not_an_event',{trial.odorName}))) = [];
        keep.laser_amps = [trial.laserAmp];
        keep.laser_pows = {trial.laserPower};
        keep.laser_times = [trial.laserTimes];
        align = 'laser';
        %both stimuli present
    case {'odor_light'}
        trial([trial.laserAmp]<1) = [];
        trial(strcmpi('none',{trial.odorName}) | strcmpi('empty',{trial.odorName}) | strcmpi('empty',{trial.odorName}) | strcmpi('not_an_event',{trial.odorName})) = [];
        keep.odors       = {trial.odorName};
        keep.concs       = [trial.odorConc];
        keep.laser_amps  = [trial.laserAmp];
        keep.laser_pows  = {trial.laserPower};
        keep.laser_times = [trial.laserTimes];
        align = 'sniff';
end

nt=numel(trial);
if nt<1
    raster = [];
    return
end

trialId = {trial.id};
tVec = (t1:t2);
spikes = zeros(nt,length(tVec));
t0Vec = nan(1,numel(trial));
x    = zeros(1e5,1);    y    = zeros(1e5,1);
nsp = 0;

for it = 1:nt
    %get all the spikes for this unit within the window tVec
    
    ii = find(trial(it).sniffZeroTimes(1,:)>trial(it).odorTimes(1)*1.009,1);
    switch align
        case 'sniff'
            if isempty(ii) || ii<1
                continue
            end
            t0 = trial(it).sniffZeroTimes(1,ii) + trial(it).start;
        case 'laser'
            t0 = trial(it).laserTimes(1) + trial(it).start;
            if isempty(ii) || ii<1
                delay = nan;
            else
                delay = trial(it).laserTimes(1) - trial(it).sniffZeroTimes(1,ii);
            end
            keep.laserDelay(it) = delay;
    end
    
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

%     figure
%     plot(x,y,'.','MarkerSize',7);
%
%%%%%%%%%%%%%%%

% keep the fields according to the stim type
raster         = structfun(@(x) x, keep, 'UniformOutput', false);
% then add everything else
raster.trialId = trialId;
raster.spikes  = spikes;
raster.t       = tVec;
raster.t0      = t0Vec;
%for quick debugging of rasters
raster.t1      = t1;
raster.t2      = t2;
raster.x = x;
raster.y = y;
end


%%% Script for creating big raster of all the trials for a cell for neil,
%%% using the trial structure from export_data
%returns one raster for one particular rec in which the cell appears.
%unitId is the unit identifier trhough the session.
function [rasters] = just_a_raster_set(mouse, sess, rec, unitSessNumber, doit, sType)
%for all the sTypes, make the rasters for this cell
%mouse
%sess
%rec
%unitSessNumber
%doit (0, 1, 2) (don't, for neil, for further analyisis in python)

%sType          : cell array of stim types. default is {'odor' 'light'
%'odor_light'}
if nargin<6 || isempty(sType)
    sTypes = {'odor' 'light' 'odor_light'};
end
sTypes = {'odor'};
%unit is a unit of the type of neil units
%get the data of the unit, get all the trials the unit is in (from neil
%trials, and make a raster centered on the first inhale after onset of
%odor
%
%unitNumber is which unit of that rec you want to get the raster.
fprintf('Making rasters for unit %s_%03d_%03d_%s\n',mouse,sess,unitSessNumber,rec)
fn = file_names(mouse,sess,rec);
trial = neil_trial_structure(mouse,sess,rec, doit);

if all(strcmpi(sTypes,{'odor'}))
    non_odor_trials = (cellfun(@(x) any(strcmpi(x, {'none', 'empty', 'not-an-event', 'dummy'})),{trial.odorName}));
    trial(non_odor_trials)=[];
end

save(fn.exp_trial,'trial');
load(fn.exp_spikes);

%meta data
load(fullfile(fn.fold_exp_data,'unitsmeta.mat'));
%look up the cell
cellsOfRec =cellsArray( find( strcmpi(mouse,{cellsArray.mouse}) & strcmpi(rec,{cellsArray.rec}) & [cellsArray.sess]==sess));
thisCell = cellsOfRec([cellsOfRec.sessCell]==unitSessNumber);
unitNumber = find(strcmpi(thisCell.uId,{cellsOfRec.uId})); %number of cell amongst the cells of the rec
thisUnit = unit(unitNumber);
%save the unit to a single file (it makes it easier to read)
if doit == 2
    unitFileName = fullfile(fn.fold_exp_data, sprintf('%s_%03d_%s_%03d_spikes.mat', mouse, sess, rec,thisUnit.sessCell));
    save(unitFileName, 'thisUnit');
end
%     %select a particular odor, for debugging purposes
%     trial=trial(strcmpi('2-hydroxyacetophenone',{trial.odorName}));
%     nt = numel(trial);


%the trials in exp are already refined (output of neil_trial_structure)
%get that trial structure and go trial by trial adding rows to the big
%raster
for i=1:numel(sTypes)
    sType = sTypes{i};
    a_raster = do_raster(trial, round(thisUnit.times),sType);
    if ~isempty(a_raster)
        a_raster.rec     = rec;
        a_raster.cell    = thisCell;
        rasters.(sType) = a_raster;
    end
end

end

function [raster] = just_a_baseline(mouse,sess,rec,unitSessNumber)
%%% Script for creating big raster of all the trials for a cell for neil,
%%% using the trial structure from export_data
%returns one raster for one particular rec in which the cell appears.
%unitId is the unit identifier trhough the session.
fprintf('Making baseline for unit %s_%03d_%03d_%s\n',mouse,sess,unitSessNumber,rec)
fn = file_names(mouse,sess,rec);

trial = assemble_baseline_trial_tructure_zk(mouse,sess,rec);
load(fn.exp_spikes);

%meta data
load(fullfile(fn.fold_exp_data,'unitsmeta.mat'));
%look up the cell
cellsOfRec =cellsArray( find( strcmpi(mouse,{cellsArray.mouse}) & strcmpi(rec,{cellsArray.rec}) & [cellsArray.sess]==sess));
thisCell = cellsOfRec([cellsOfRec.sessCell]==unitSessNumber);
unitNumber = find(strcmpi(thisCell.uId,{cellsOfRec.uId})); %number of cell amongst the cells of the rec
thisUnit = unit(unitNumber);


%the trials in exp are already refined (output of neil_trial_structure)
%get that trial structure and go trial by trial adding rows to the big
%raster
t1 = 1;
t2 = 3000;

nt=numel(trial);
if nt<1
    raster = [];
    return
end

trialId = {trial.trialUId};
tVec = (t1:t2);
spikes = zeros(nt,length(tVec));
t0Vec = nan(1,numel(trial));
sp = round(thisUnit.times);

x    = zeros(1e5,1);    y    = zeros(1e5,1);
nsp = 0;

for it = 1:nt
    t0 = trial(it).start;
    
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

%     figure
%     plot(x,y,'.','MarkerSize',7);
%
%%%%%%%%%%%%%%%
% Add cellId?  mouse_sess_rec_int --> int=unique integer of cell
%                                         in session from unitDb
raster.cellId  = thisCell.Id;
raster.spikes  = spikes;
raster.trialId = trialId;
raster.mouse   = mouse;
raster.sess    = sess;
raster.rec     = rec;
raster.uid     = thisCell.uId;
raster.t0      = t0Vec;
raster.t1      = t1;
raster.t2      = t2;
%for quick debugging of rasters
raster.x = x;
raster.y = y;

end

function trials = neil_trial_structure(mouse,sess,rec, doit)
%gets a trial structure for a mouse,sess,rec and writes it in a format
%suitable for neil
fn=file_names(mouse,sess,rec);
load(fn.trial);
load(fn.sess_info);

recInfo =  info.rec(strcmpi(rec,{info.rec.name}));
rsmSniff = load(fn.rsm_data,'Sniff');

%pick all the trials and get odor, laser and sniff info.
%this will be used to produce the rasters for light, odor and light+odor
%trials

t1 = -200;
t2 = 2500-1;

trials = [];
for it=1:numel(trial)
    ttr = trial(it);
    tr.start      = ttr.start+t1;
    if isempty(tr.start) || isnan(ttr.start) || ttr.start+t1<1
        continue
    end
    tVec          = ttr.start + (t1:t2-1); %absolute times of the trial chunk
    tr.id         = sprintf('%s_%03d_%s_%d',mouse,sess,rec,ttr.start);
    tr.odorName   = ttr.odorName;
    tr.odorConc   = ttr.odorConc;
    %laser data
    tr.laserAmp   = ttr.laserAmp;
    tr.laserPower = ttr.laserPower;
    tr.laserTimes = ttr.laserTimes-t1;
    
    %the cutout of the sniff trace
    try
    tr.flow       = rsmSniff.Sniff(tVec);
    catch me
        warning('bad sniff extract')
    end
    %the stimulus vector
    tr.stim       = zeros(1,t2-t1);
    tr.odorTimes  = ttr.odorTimes - t1 ;
    stimTimes  = tr.odorTimes;
    stimTimes([tr.odorTimes]>t2-t1)=t2-t1;
    stimTimes([tr.odorTimes]<0)=0;
    if all(stimTimes>0)
        tr.stim(stimTimes(1):stimTimes(2))=1;
    end
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
            spZeros(2,:) = round((ttr.sniffParabZeroTimes(2,:)))-t1;
            spZeros(1,:) = round((ttr.sniffZeroTimes(1,:)))-t1;
            tr.sniffZeroTimes = ttr.sniffZeroTimes - t1;
            if doit ==2
                tr.spZeros = spZeros;
            end
            %         sz = spZeros - t1*ones(size(spZeros));
            %         firstSnif = find(sz>0,1,'first');
            %         lastSnif  = find(sz<t2,1,'last');
            %         sn(1:sz(firstSnif))=bitget(firstSnif-1,1);
            %         for is=firstSnif:lastSnif-1
            %             sn(sz(is):sz(is+1))=bitget(is,1);
            %         end
            %         sn(sz(lastSnif):end)=bitget(lastSnif,1);
            %         sn(sn==0)=-1;
            %         tr.sniffPhase = sn;
            trials = [trials tr];
        end
    catch
        warning('Error in sniff phases of trial %d - %d',ttr.run,ttr.runTrialNum);
        tr.sniffPhase = nan;
    end
end

end

%Kristina's functions to do the baselines:
function [trialsBase] = assemble_baseline_trial_tructure_zk(mouse,sess,rec)
fn=file_names(mouse,sess,rec);
% load trial structure
q = load(fn.trial);
trial=q.trial; clear q;

% load sniffs
q = load(fn.sniffs);
sniff = q.sniff;

q = load(fn.rsm_data);
SniffTrace = q.Sniff;

% get Fvpin on
%fv.name   = 'finalValve';
%fv.type   = 'digital';
%fv.chanId = 'fvPin';
%fv.events = get_analog_events(fv.chanId,mouse,sess,rec,run,'figures','noplot');

%%clean up trials for sessions < 20 (before serial trial number was used)
badTrials=arrayfun(@(x) isempty(trial(x).start),1:numel(trial));
if any(badTrials)
    warning('Some bad trials found in the rec (start was empty)');
end
trial(badTrials)=[];

%FVpin = q.FVpin;
%clear q;
%finalValves on/off times:
fvOnTimes = [[trial.odorTimes]' [trial.start]'];
fvOnTimes(fvOnTimes(:,1)==0,:)=[];
fvOnTimes(:,1)=fvOnTimes(:,3)+fvOnTimes(:,1);
fvOnTimes(:,2)=fvOnTimes(:,3)+fvOnTimes(:,2);

laserOnTimes = [[trial.laserTimes]' [trial.start]'];
laserOnTimes(laserOnTimes(:,1)==0,:)=[];
laserOnTimes(:,1)=laserOnTimes(:,3)+laserOnTimes(:,1);
laserOnTimes(:,2)=laserOnTimes(:,3)+laserOnTimes(:,2);

trialsBase = [];
badSniffs = 0;
%get all the sniffs that are before an open valve
eventOnTimes = [fvOnTimes(:,1)];
for ifv=1:numel(eventOnTimes)
    %for now, forget about the sniff zeros and such
    %just use the time stamps of the fv opening to chop the sniff
    
    %find the sinff before a trial, light or odor
%     prev_sniff = find([sniff.t0]< eventOnTimes(ifv,1)-3000,1,'last');
%     if(isempty(prev_sniff))
%         continue
%     end
%     sn=sniff(prev_sniff);
    %t_inh = sn.t0+sn.t_zer(1);
    
    t_inh = eventOnTimes(ifv,1);
    t1= -3000;
    t2= 0;
    
    % is a sniff not within a stimulus
    %sn_zeros = sn.t_zer - (150+t1);
    
    tr.start = t_inh +t1;
    
    %if any(sn_zeros(2:3)<1) || (t_inh + t2) > numel(SniffTrace)
    if (t_inh + t2) > numel(SniffTrace)
        %warning('problem with sniff zeros at sn.t0 %d',sn.t0)
        badSniffs = badSniffs+1;
        continue
    end
    tr.trialUId = [fn.basename_an 'trial' num2str(tr.start)];
    try
        tr.sniffFlow = SniffTrace(tr.start:t_inh+t2);
    catch me
        warning('sniff begins too early or ends too late for chopping the right segment around it (sniff %d)',ifv);
        continue;
    end
    %tr.sniffPhase = ones(size(tr.sniffFlow));
    %tr.sniffPhase(1:-t1)=-1;
    %tr.sniffPhase(sn_zeros(2):sn_zeros(3)) = -1;
    trialsBase = [trialsBase tr];
    
    %if it got sniff, 
    
    %for debugging
    %         figure
    %         plot(-SniffTrace(tr.start:(tr.start+1000)));
    %         hold on
    %         plot(-tr.sniffFlow,'r');
    %         plot(sn_zeros,[0 0 0],'b*');
    %         plot(tr.sniffPhase*100);
    
end
if badSniffs>0
    warning('%d problems with sniff zeros',badSniffs);
end
fprintf('(%d sniffs)\n', numel(trialsBase));
save(fullfile(fn.fold_exp_data,sprintf('%strialsBase.mat',fn.basename_an)), 'trialsBase');

end

function [trialsBase] = assemble_baseline_trial_tructure(mouse,sess,rec)
fn=file_names(mouse,sess,rec);
% load trial structure
q = load(fn.trial);
trial=q.trial; clear q

% load fields from rsm file
q = load(fn.rsm_data);
Sniff = q.Sniff;
FVpin = q.FVpin;
clear q

trialsBase=struct;
otr=0;
for itr=1:numel(trial)
    tr=trial(itr);
    
    if ~isempty(tr.start) && ~isempty(tr.sniffParabZeroTimes)
        if size(tr.sniffParabZeroTimes,2) == size(tr.sniffZeroTimes,2)
            sniffZeros_toUse = nan(2,size(tr.sniffParabZeroTimes,2));
            sniffZeros_toUse(2,:) = round(tr.sniffParabZeroTimes(2,:));
            sniffZeros_toUse(1,:) = tr.sniffZeroTimes(1,:);
        else
            disp(itr)
            error('sniff zeros and parabs are not same length')
        end
        
        % make trial unique id
        trialUId = [fn.basename_an 'trial' num2str(tr.start)];
        
        if diff([tr.odorTimes])>0 && tr.odorConc>0
            
            % get 3000ms of sniff data pre FV onset
            FV_recTime = round(tr.odorTimes(1)*1.009) + tr.start;
            
            t2 = FV_recTime;
            t1 = t2-3000;
            if t1<1
                continue
            end
            %         itr
            %         if itr~=1 && t1<(trial(itr-1).start+trial(itr-1).runTrialDur)
            %             t1=t2-1000;
            %         end
            
            FV_rsmCandidates =  max(1,tr.start-2e4) + find(diff(FVpin(max(1,tr.start-2e4):min(tr.start+2e4,end)))>1e4);
            %                 FV_rsmCandidates =  (max(1,tr.start-2e4) + find(diff(FVpin(max(1,tr.start-2e4):tr.start+2e4))>1e4)+1);
            [~,FVtime_idx] =  min(abs(FV_rsmCandidates-FV_recTime));
            FV_rsmTime = FV_rsmCandidates(FVtime_idx);
            
            ms_tr2rsm = diff([FV_recTime FV_rsmTime]);
            if abs(ms_tr2rsm)>3
                %                 error('Couldnt find close enough FVon time in rsm_data')
                continue
            else
                otr=otr+1;
                % Create sniff phase vector
                vectorizedSniff = reshape(sniffZeros_toUse+tr.start,1,2*size(sniffZeros_toUse,2));
                vectorizedSniff = vectorizedSniff(vectorizedSniff>t1&vectorizedSniff<t2)-t1;
                
                inhOnsets = sniffZeros_toUse(1,:)+tr.start > t1 & sniffZeros_toUse(1,:)+tr.start < t2;
                pauseOnsets = sniffZeros_toUse(2,:)+tr.start > t1 & sniffZeros_toUse(2,:)+tr.start < t2;
                
                phase_inRange = [inhOnsets; -1.*pauseOnsets];
                vectorizedPhase = phase_inRange(phase_inRange~=0)';
                
                if size(Sniff(t1+1+ms_tr2rsm:t2+ms_tr2rsm),2) ~=3000
                    aaa=243;
                end
                
                if numel(vectorizedSniff)<5
                    warning('Skipped a trial with few sniffs detected')
                    trialsBase(otr).trialId = trialUId;
                    trialsBase(otr).start = t1;
                    trialsBase(otr).sniffFlow = Sniff(t1+1+ms_tr2rsm:t2+ms_tr2rsm);
                    trialsBase(otr).sniffPhase = nan;
                else
                    
                    if size(vectorizedSniff)==size(vectorizedPhase)
                        vectorizedPhase = [vectorizedPhase(2) vectorizedPhase];
                        vectorizedSniff = [1 vectorizedSniff t2-t1];
                    else
                        error('sniff vectors not same size')
                    end
                    if ~any(diff(vectorizedSniff)>1200) && vectorizedSniff(1)<1200
                        sniffPhase = nan(1,t2-t1);
                        for ib = 1:length(vectorizedSniff)-1
                            sniffPhase(vectorizedSniff(ib):vectorizedSniff(ib+1)) = vectorizedPhase(ib);
                        end
                        
                        ii = find(tr.sniffZeroTimes(1,:)>tr.odorTimes(1)*1.009,1);
                        figure(1); clf; hold on
                        plot(1:t2-t1,Sniff(t1+ms_tr2rsm:t2+ms_tr2rsm-1))
                        plot(vectorizedSniff,zeros(size(vectorizedSniff)),'r*')
                        %         fill([t1;t2]*ones(1,2),ones(2,1)*[-3e4 3e4],'-','LineWidth',100,'Color',[1 0.9 0.95])
                        %                         fill([t1 t1 t2 t2],[-3e4 3e4 3e4 -3e4],[1 0.95 0.97],'EdgeColor','none')
                        %                         plot(t1+ms_tr2rsm:t2+ms_tr2rsm+1e3,Sniff(t1+ms_tr2rsm:t2+ms_tr2rsm+1e3))
                        plot((tr.sniffZeroTimes(1,ii)+tr.start)*ones(1,2),[-3e4 3e4],':g','LineWidth',3)
                        plot((t1+1:t2),-1e4.*sniffPhase,'k')
                        pause(2)
                        
                        trialsBase(otr).trialId = trialUId;
                        trialsBase(otr).start = t1;
                        trialsBase(otr).sniffFlow = Sniff(t1+1+ms_tr2rsm:t2+ms_tr2rsm);
                        trialsBase(otr).sniffPhase = sniffPhase;
                    else
                        trialsBase(otr).trialId = trialUId;
                        trialsBase(otr).start = t1;
                        trialsBase(otr).sniffFlow = Sniff(t1+1+ms_tr2rsm:t2+ms_tr2rsm);
                        trialsBase(otr).sniffPhase = nan;
                        
                        %                         figure(1); clf; hold on
                        %                         plot(1:t2-t1,Sniff(t1+ms_tr2rsm:t2+ms_tr2rsm-1))
                        %                         plot(vectorizedSniff,zeros(size(vectorizedSniff)),'r*')
                        warning('Sniff_analysis may have missed a cycle; check waveform')
                    end
                end
            end
        end
    end
end
% save(fullfile(fn.base_folder,'data_Neil',sprintf('%strialsBase.mat',fn.basename_an)), 'trialsBase')
end

function [spikesBase,qCells] = assemble_baseline(qCells)
% This function creates a file of cell/spike data formatted for Neil.
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
    
    % Extract spikes and save structure
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
    save(fullfile(fn.base_folder,'data_Neil',sprintf('%s%i_spikesBase.mat',fn.basename_an,qCells(ic).sessCell)), 'spikesBase')
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

function [noStimSniffs] = no_stim_sniffs_raster(unit_meta)
% get the sniffs with no stimulus.
% save the sniff traces, the inhales and exhales ref. to the beginning of
% the rec.

cp = cell_passport_tools();
noStimSniffs = cp.get_no_stim_sniffs(unit_meta);

end



