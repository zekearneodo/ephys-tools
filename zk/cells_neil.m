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

% go to the trial structure and make a copy with the adecquate spikeTimes

%get the trial structures, empty the spikes and write it in the export
%folder
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
            fprintf('   rec %s ...',rec)
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


    