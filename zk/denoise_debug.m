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
raw_file = fopen(fn.ss_rec);
ds=data_structure_tools;
%make the matrix with the shank sites

%get the info
q=load(fn.ss_sess_info);
[rec_info,~]=ds.get_info(mouse,sess,rec,1);

chan_group='d';
group_info = rec_info.chGroupInfo(strcmpi(chan_group,[rec_info.chGroupInfo.name]))
chan_list = group_info.chan+1;

%load the channels into a matrix
rf=dir(fn.ss_rec);
data_points = rf.bytes/2/rec_info.nChan;
raw_data = zeros(data_points,group_info.nCh);
for i=1:numel(chan_list)
    ch = chan_list(i);
    fprintf('Reading site %d (total %d sites)\n',ch,numel(chan_list));
    raw_data(:,i) = pp.read_analog_channel(raw_file,chan_list(i),rec_info.nChan);
end


%plot a segment
chunk = 1:100000;
data_range=(chunk);
plot_range = [1:19530];
data = raw_data(data_range,:);
data=bsxfun(@minus,data,mean(data));

fig=pp.plot_neural_data(data(plot_range,:));

%start playing with the segment
[u,s,v]=svd(data);


%substract all the common components to one channel
data_noise = zeros(size(data));
data_clean = data;
for ich=1:numel(chan_list);
    chans=1:numel(chan_list);
    chans_noti=chans(~(chans==ich));
    %pp.plot_neural_data(data(plot_range,chans_noti));
    [u,s,v]=svd(data(:,chans_noti),'econ');
    b=zeros(1,3);
    for ic=1:3
        b(ic)=dot(data(:,ich),u(:,ic))/norm(u(:,ic));
        data_noise(:,ich)=data_noise(:,ich)+b(ic)*u(:,ic);
    end
    data_clean(:,ich)=data_clean(:,ich)-data_noise(:,ich);
end
pp.plot_neural_data(data_noise(plot_range,:));
pp.plot_neural_data(data_clean(plot_range,:));

% figure
% plot(data(:,ich))
% hold
% plot(data_clean(:,ich))

%get the spike candidates and replace the supra-trheshold events by the
%common noise component
%get spikes in the cleaned matrix
spikes = pp.detect_spikes(data_clean,p);
%replace spikes with noise data in a copy of the raw data
data_nospikes=data;
for ich = 1:numel(chan_list)
    chans = 1:numel(chan_list);
    chans_noti = chans(~(chans==ich));
    events     = spikes{ich};
    ev_select  = cell2mat(arrayfun(@(x) (x-p.extract_s_before):(x+p.extract_s_before),events,'UniformOutput',false));
    data_nospikes(ev_select,ich) = data_noise(ev_select,ich);
end

% now do the second stage:
% go through the data and subtract the projections of the pcs
% but now the data_noti comes from data_nospikes

%substract all the common components to one channel
data_noise = zeros(size(data));
data_clean = data;
for ich=1:numel(chan_list);
    chans=1:numel(chan_list);
    chans_noti=chans(~(chans==ich));
    %pp.plot_neural_data(data(plot_range,chans_noti));
    [u,s,v]=svd(data_nospikes(:,chans_noti),'econ');
    b=zeros(1,3);
    for ic=1:3
        b(ic)=dot(data(:,ich),u(:,ic))/norm(u(:,ic));
        data_noise(:,ich)=data_noise(:,ich)+b(ic)*u(:,ic);
    end
    data_clean(:,ich)=data_clean(:,ich)-data_noise(:,ich);
end
pp.plot_neural_data(data_noise(plot_range,:));
pp.plot_neural_data(data_clean(plot_range,:));
