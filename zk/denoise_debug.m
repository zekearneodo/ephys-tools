% Program for testing common pca removing before spike sorting
%
global pp;
post_exp_processing_voyeur_ndm_par
%
mouse = 'ZKawakeM72';
sess  = 20;
rec   = 'c';
fn=file_names(mouse,sess,rec);

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
data_range=(1:1953);
data = raw_data(data_range,:);

fig=pp.plot_neural_data(data)

%start playing with the segment
[u,s,v]=svd(data);
pp.plot_neural_data(u(:,1:8));
pp.plot_neural_data(u(1:8,:)');

    data_clean=data;
    data_noise=zeros(size(data));
%substract all the common components to one channel
for ich=1:8;
    chans=1:numel(chan_list);
    chans_noti=chans(~(chans==ich));
    data_noti=data(:,chans_noti);

    %pp.plot_neural_data(data_noti);
    [u,s,v]=svd(data_noti,0);
    b=zeros(1,3);
    for ic=1:3
        b(ic)=dot(data(:,ich),u(:,ic))/norm(u(:,ic));
        data_noise(:,ich)=data_noise(:,ich)+b(ic)*u(:,ic);
    end
    data_clean(:,ich)=data_clean(:,ich)-data_noise(:,ich);
end
pp.plot_neural_data(data_noise);
pp.plot_neural_data(data_clean);

figure
plot(data(:,ich))
hold
plot(data_clean(:,ich))

