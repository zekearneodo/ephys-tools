function fill_x_run(mouse, sess,rec,runInfo,group,data)
%write an array of vector variables as datasets in the 'type' group
%of the corresponding hdf5 file into the corresponding rec, run.
%input:
%mouse
%sess
%rec
%runInfo : info structure for the run
%group   : type of the dataset ('streams','events')
%data    : structure with vectors to append as datasets to the group
%    .name_of_dataset : dataset

fn=file_names(mouse,sess,rec);
q=load(fn.sess_info);
sInfo=q.info;

fi=fields(data);
fieldsToWrite=fi(cellfun(@(x) ~any(strfind(x,'meta')),fi));
%
parentGroup=sprintf('/%s/%02d/%s',rec,runInfo.num,group);
for iF=1:numel(fieldsToWrite)
    %make a dataset with the name and size and write the data
    field=fieldsToWrite{iF};
    nPoints=numel(data.(field));
    chunkSize = min(1024,nPoints);
    h5create(fn.sess_h5,parentGroup,[1 nPoints],'Datatype',class(data.(field)),'ChunkSize',[1 chunkSize]);
    h5write(fn.sess_h5,parentGroup,data.(field),[1 1],[1 nPoints]);
end

end
