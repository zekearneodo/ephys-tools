%very basic script to remove the first candidates of common noise
%step 0:
%find the spike candidates that appear exactly at the same time in all the
%groups, and make them a new cluster.
%step 1:
%instead of removing them from the pool, substract the average in all the
%ephys channels; then run the detection again
%step 2:
%get the raw data file, find the PCs for the groups at the spike
%candidates. If the first PC is the same for all the groups, substract the
%projection to each individual channel.

%Begin with step 0:
%Open the fet files:
mouse='SDbehavingM72';
sess=1;
rec='a';
irec=1;
jitter=(-10:10);

fn=file_names(mouse,sess,rec);
q=load(fn.ss_sess_info);
thisRec=q.info.rec(irec);
%list and load the .res files (the timestamps of spike candidates)
resFiles=dir(fullfile(fn.fold_ss_rec,['rec_' rec '.res.*']));
cluFiles=dir(fullfile(fn.fold_ss_rec,['rec_' rec '.clu.*']));

for irf=1:numel(resFiles)
    group(irf).res=load(fullfile(fn.fold_ss_rec,resFiles(irf).name),'-ascii');
    group(irf).len=numel(group(irf).res);
    group(irf).clu=load(fullfile(fn.fold_ss_rec,'CLU_00',cluFiles(irf).name),'-ascii');
    group(irf).nClu=group(irf).clu(1);
end

[~,lenSorted]=sort([group.len]);

%find recursively the elements of the shorter that appear in the longer
%(within 2 samples)
common=group(lenSorted(1)).res;

for i=2:numel(group)
    ig=lenSorted(i);
    lCommon=zeros(size(common));
    for jit=jitter       
        lCommon=lCommon+ismember(common,group(ig).res+jit);
    end
    common=common(lCommon>0);
end


%update the cluster file to:
%have one more cluster
%have all the common belong to the new cluster
for ig=1:numel(group)
    
    lthisGroup=zeros(size(group(ig).res));
    %the (logical) indexes of the spikes of this group that appear in the found
    %common spikes
    for jit=jitter       
        lthisGroup=lthisGroup+ismember(group(ig).res+jit,common);
    end
    
    
    group(ig).clu(1)=group(ig).nClu+1;
    group(ig).clu(find(lthisGroup)+1)=group(ig).nClu+1;
    fId=fopen(fullfile(fn.fold_ss_rec,[cluFiles(ig).name]),'w');
    fprintf(fId,'%d\n',group(ig).clu);
    fclose(fId);
end

