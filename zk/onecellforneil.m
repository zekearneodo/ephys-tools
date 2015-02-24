% one unit for neil.
% - calls the function just_a raster, that makes the trial structure and the
%   spikes matrix
% - makes the raster from the spikes matrix to compare with the rasters we
%   already have and check consistency (i.e.; timestamps not screwed up).
% 

mouse    = 'ZKawakeM72';
sess     = 13;
rec      = 'f';
sessCell = 10;

cn=cells_neil();
fn=file_names(mouse,sess,rec);
% trial = neil_trial_structure(mouse,sess,rec);
% save(fn.exp_trial);
raster = cn.just_a_raster(mouse,sess,rec,sessCell);


% plot the raster for an odor
% get the rasters from the spikes, to confirm that the spikes matrix is ok

odor = '2-hydroxyacetophenone';
odorTrials=find(strcmpi(odor,raster.odors));
odorSpikes = raster.spikes(odorTrials,:);

nt=numel(odorTrials)
x=[];
y=[];
t1= -200;
for it = 1:nt
    %disp(it)
    spikes = find(odorSpikes(it,:));
    n = numel(spikes);
    x = [x spikes];
    y = [y it*ones(1,length(spikes))];
end
figure
plot(x,y, '.', 'MarkerSize',7)
%figure
%plot(raster.x,raster.y, 'r.', 'MarkerSize',7)

%%% make the cell structure
%the raster

%the cell data
% get the sniff and check the t0
rsmSniff = load(fn.rsm_data,'Sniff');
sniff=rsmSniff.Sniff*(-1.);
figure
hold
plot(sniff)
t0=raster.t0;
plot(t0,zeros(size(t0)),'r*')

%get all segments 200 ms before and 200 ms after and plot them all together
ts1 =-200;
ts2 = 200;

sniff_cuts=zeros(numel(t0),(ts2-ts1+1));
for i=1:numel(t0)
    sniff_cuts(i,:) = sniff(t0(i)+ts1:t0(i)+ts2);
end
figure
plot(sniff_cuts')