%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script for creating big raster of all the trials for a cell for neil,
%%% using the trial structure from export_data

function [raster] = raster_neil(mouse,sess,rec,unitNumber)
    %unit is a unit of the type of neil units
    %get the data of the unit, get all the trials the unit is in (from neil
    %trials, and make a raster centered on the first inhale after onset of
    %odor
    % 
    %unitNumber is which unit of that rec you want to get the raster.
    
    fn = file_names(mouse,sess,rec);
    q=load(fn.exp_spikes);
    unit = q.unit;
    q=load(fn.exp_trial);
    trial = q.trial;
    
%     %select a particular odor, for debugging purposes
%     trial=trial(strcmpi('2-hydroxyacetophenone',{trial.odorName}));
%     nt = numel(trial);

    
    %the trials in exp are already refined (output of neil_trial_structure)
    %get that trial structure and go trial by trial adding rows to the big
    %raster
    t1 = -200;
    t2 = 400;
    
    nt=numel(trial);
    odors   = {trial.odorName};
    concs   = [trial.odorConc];
    trialId = {trial.id};
    tVec = (t1:t2);
    spikes = zeros(nt,length(tVec));
    sp = round(unit(unitNumber).times);
    
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
        %
        if it==95
          disp(it);
        end
        disp(it);
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
    raster.odors = odors;
    raster.concs    = concs;
    raster.trialId = trialId;
    raster.spikes  = spikes;
    raster.t       = tVec;
    raster.rec     = rec;
    raster.unitNumber = unitNumber;
    
    %for quick debugging of rasters
    raster.x = x;
    raster.y = y;
    
    %data for the unit
    cell.mouse  = mouse;
    cell.sess   = sess;

    %meta data
    load(fullfile(fn.fold_exp_data,'unitsmeta.mat'));
    %look up the cell
    cellsOfRec = find( strcmpi(mouse,{cellsArray.mouse}) & strcmpi(rec,{cellsArray.rec}) & [cellsArray.sess]==sess);
    thisCell = cellsArray(cellsOfRec(unitNumber))
    cell.uid = sprintf('%s_%03d_%03d',mouse,sess,thisCell.sessCell);
    cell.light = thisCell.light;

    cell.raster = raster;
    
    %save it
    cellFn=fullfile(fn.fold_exp_data,sprintf('%s_cell.mat',cell.uid));
    save(cellFn,'-struct','cell');
end