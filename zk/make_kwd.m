function outFile = make_kwd(mouse,sess,rec)
%but no; I am not going to generate this fle.
%I will pass the 

% create the file for the rec, the group for the run
% set the attributes
% read one channel and put it in the corresponding dataset.
[fn,recInfo] = file_names(mouse,sess,rec);
commonAttributes = struct('CLASS','GROUP','TITLE','','VERSION','1.0');

nChans  = 5;
nPoints = 20000*120;

%remove the file
if exist(h5fp)
   delete(h5fp)
end
h5create(h5fp,'/recordings/0/data',[nChans nPoints],'Datatype','int16','ChunkSize',[nChans 1024])
h5writeatt(h5fp,'/recordings/0/data','VERSION',1.1)
h5writeatt(h5fp,'/recordings/0','downsample_factor','N.')
%h5writeatt(h5fp,'/recordings/0/filter','CLASS','GROUP')

%put one channel into the dataset
chan = 1;
%read the channel 
runRawFile = fullfile(fn.fold_rd_sess,sprintf('%s_%02d_data.bin',rec,run));
fid=fopen(runRawFile,'r');
fseek(fid, 2*(chan-1), 'bof');
Y = fread(fid, inf, 'int16', (nChans-1)*2);
fclose(fid);
%write it in the dataset
Y = int16(transpose(Y(1:nPoints)));
h5write(h5fp,'/recordings/0/data',Y,[1 1],[1 nPoints])

    function make_file_structure
    %make the empty fie with the containers (groups, datasets)
    %Make the file
    fcpl = H5P.create('H5P_FILE_CREATE');
    fapl = H5P.create('H5P_FILE_ACCESS');
    h5fid = H5F.create(fn.ss_raw_kwd,'H5F_ACC_TRUNC',fcpl,fapl);
    %Make the group structure
    
    
    H5F.close(h5fid);
    end
end