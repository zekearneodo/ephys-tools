function cell = merge_clusters(mouse,sess,rec,cluList)
% merges a list of cluster structures into one single cell structure
% concatenates and sorts all the timestamps and list of clusters,
% copies the features that are common to all.
fn = file_names(mouse,sess,rec);
load(fn.spikes);

%get the units from that cluster
%unitsIdx=arrayfun(@(x) find([unit.clu]==x,1),cluList);
units=unit(cluList);


if numel(unique([units.group]))>1
    error('Cant treat clusters of different channel groups as the same unit')
end

%copy the fields to be kept
keepFields = {'group','nChans','chans','sites'};
for iF = 1:numel(keepFields)
        cell.(keepFields{iF})=unit(1).(keepFields{iF});
end

%merge the fields to be merged
catFields = {'pkStamps','times','clu'};
for iF=1:numel(catFields)
    f=catFields{iF};
    fContent=arrayfun(@(x) transpose(x.(f)),units,'UniformOutput',false);
    cell.(f)=sort(cell2mat(fContent));
end

end