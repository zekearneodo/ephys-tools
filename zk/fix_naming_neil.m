%quick fix for naming mismatch between penikis baseline files and mine
% KP's files were one per cellId (i.e, one per cell, per rec)
% Gather them by uid, and save a file per cell, with an array of rasters
% (one per rec the cell appears in)
% ZK - 03-18-2015

function fix_naming_neil()
fn=file_names();
cells_folder=fullfile(fn.fold_exp_data,'neil','all_cells');

base_files = dir(fullfile(cells_folder,'*spikesBase.mat'));
all_cells=[];
for i=1:numel(base_files)
    q=load(fullfile(cells_folder,base_files(i).name));
    id_parts=strsplit(q.spikesBase.cellId,'_');
    q.spikesBase.mouse = id_parts{1};
    q.spikesBase.sess  = str2num(id_parts{2});
    q.spikesBase.rec   = id_parts{3};
    cell_n = str2num(id_parts{4});
    q.spikesBase.uid = sprintf('%s_%03d_%03d',q.spikesBase.mouse,q.spikesBase.sess,cell_n);
    all_cells=[all_cells q.spikesBase];
end

uids=unique({all_cells.uid});
for iu=1:numel(uids)
    %get all the units with the same uid
    uid = uids{iu};
    instances = all_cells(strcmpi({all_cells.uid},uid));
    disp(numel(instances))
    %make one structure with the spikes of all the instances
    cell_base.mouse = instances(1).mouse;
    cell_base.sess  = instances(1).sess;
    cell_base.uid   = uid;
    cell_base.raster = [];
    for ii=1:numel(instances);
        instance = instances(ii);
        base_raster.trialId = instance.trialId;
        base_raster.spikes  = instance.spikes;
        base_raster.rec = instance.rec;
        cell_base.raster = [cell_base.raster base_raster];
    end
    %load the cell and append the base raster
    name = sprintf('%s_spikesBase.mat',instances(1).uid);
    spikesBase=instances
    save(fullfile(cells_folder,name),'spikesBase')
end

%now remove all the bad named spikesBase files
for i=1:numel(base_files)
    delete(fullfile(cells_folder,base_files(i).name));
end
end