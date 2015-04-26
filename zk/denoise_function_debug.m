% Program for testing common pca removing before spike sorting
%
global pp;
post_exp_processing_voyeur_ndm_par
%
mouse = 'ZKawakeM72';
sess  = 20;
rec   = 'c';
fn=file_names(mouse,sess,rec);

%parameters
sample_rate = 19385;
p.nexcerpts = 50;
p.excerpt_size = round(1*sample_rate); %1 second
p.threshold_strong_std_factor = 4;
p.threshold_weak_std_factor = 2;
p.detect_spikes = 'negative';
p.extract_s_before = 16;
p.extract_s_after  = 16;

rate  = 19385;
low   = 0.100;
high = 19385/2.1;

[filt_a,filt_b] = butter(3,low*2/rate,'high');


%open the file

ds=data_structure_tools;
%make the matrix with the shank sites

%get the info
q=load(fn.ss_sess_info);
[rec_info,~]=ds.get_info(mouse,sess,rec,1);

chan_groups={'b','c','d'};
for ig=1:numel(chan_groups)
    chan_group = chan_groups{ig}
group_info = rec_info.chGroupInfo(strcmpi(chan_group,[rec_info.chGroupInfo.name]))
chan_list = group_info.chan+1;

%load the channels into a matrix
rf=dir(fn.ss_rec);
data_points = rf.bytes/2/rec_info.nChan;
raw_data = zeros(data_points,group_info.nCh);
raw_file = fopen(fn.ss_rec);
for i=1:numel(chan_list)
    ch = chan_list(i);
    fprintf('Reading site %d (total %d sites)\n',ch,numel(chan_list));
    raw_data(:,i) = pp.read_analog_channel(raw_file,chan_list(i),rec_info.nChan);
end
fclose(raw_file);

%go by chunks calling the function
chunk_size = 100000;
chunk_overlap = round(19385/2);
n_chunks = ceil(data_points/chunk_size);
data = zeros(size(raw_data));

fprintf('Removing common pca components: %d chunks:\n',n_chunks);
fprintf('Chunk:');
for i=1:n_chunks
    if mod(i,20)==1
        fprintf(' %d',i);
    end
    if i==1
        ch_overlap = 0;
    else
        ch_overlap = chunk_overlap;
    end
    if i==n_chunks
        ch_size = data_points - (n_chunks-1) * chunk_size;
    else
        ch_size = chunk_size;
    end
    ch_start = (i-1)*chunk_size+1;
    ch_end   = ch_start + ch_size-1;
    data_temp = pp.pca_denoise_chunk(raw_data((ch_start-ch_overlap):ch_end,:),p);
    data_temp(1:ch_overlap,:)=[];
    data(ch_start:ch_end,:) =  data_temp;
end
fprintf('\n');

%save the data in the corresponding channels
sess  = 920;
rec   = 'c';
fn=file_names(mouse,sess,rec);
fid=fopen(fn.ss_rec,'r+');
for i=1:numel(chan_list)
    ch = chan_list(i);
    fprintf('Writing site %d (total %d sites)\n',ch,numel(chan_list));
    pp.write_analog_channel(data(:,i), fid, 0, ch, rec_info.nChan);
end
fclose(fid);
end
