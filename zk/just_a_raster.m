%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script for creating big raster of all the trials for a cell for neil,
%%% using the trial structure from export_data
%returns one raster for one particular rec in which the cell appears.
%unitId is the unit identifier trhough the session.

function [raster] = just_a_raster(mouse,sess,rec,unitSessNumber)
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
    raster.rec     = rec;
    raster.cell    = thisCell;
    
    %for quick debugging of rasters
    raster.x = x;
    raster.y = y;
    
end