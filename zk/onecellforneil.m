% one unit for neil.

mouse    = 'KPawakeM72';
sess     = 16;
rec      = 'a';
sessCell = 16;

fn=file_names(mouse,sess,rec);
trial = neil_trial_structure(mouse,sess,rec);
save(fn.exp_trial);
raster = just_a_raster(mouse,sess,rec,sessCell);


%plot the raster for an odor
odor = '2-hydroxyacetophenone';
odorTrials=find(strcmpi(odor,raster.odors))
odorSpikes = raster.spikes(odorTrials,:);

nt=numel(odorTrials)
x=[];
y=[];
t1= -200;
for it = 1:nt
    %disp(it)
    spikes = find(odorSpikes(it,:));
    n = numel(spikes)
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