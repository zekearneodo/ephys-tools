%Scanning trhough cells meta
ffn=file_names();
rasters_folder=fullfile(ffn.fold_an_data,'concentration','rasters');
% load all cells in a cell structure
%mice = {'ZKawakeM72','KPawakeM72'};
mice = {'ZKawakeM72'};
%load all the cells into an array of cells
disp(wrap_message('Getting cells and trial structures for Neil','*'));
cellsArray = [];
for is=1:numel(mice)
    cs=mice{is};
    cellFiles=dir(fullfile(ffn.fold_unit_db, [cs '*.mat']));
    allCells = arrayfun(@(x) load(fullfile(ffn.fold_unit_db,x.name)),cellFiles,'UniformOutput',false);
    cellsArray = [cellsArray [allCells{:}]];
end

%get the range
mouse = mice{1};
sess_min = 20;
sess_max = 100;
quality = 1;
light = true;
odor = true;
cellsArray([cellsArray.sess]<sess_min)=[];
cellsArray([cellsArray.sess]>sess_max)=[];
cellsArray(~([cellsArray.quality]==quality))=[];

if light
    cellsArray([cellsArray.light]==0)=[];
end

if odor
    cellsArray([cellsArray.odor]==0)=[];
end

%make the responses for all the cells
for ic=1:numel(cellsArray)
    theCell = cellsArray(ic);
    % make plot the responses
    vr=visualize_responses_new(theCell.mouse,theCell.sess,theCell.rec,theCell.clu,'odor');
    fn = file_names(theCell.mouse,theCell.sess,theCell.rec);
    resBaseName = sprintf('%s_%03d_%s_odor_units%02d_resp',theCell.mouse,theCell.sess,theCell.rec,theCell.clu);
    load(fullfile(fn.fold_an_sess,[resBaseName '.mat']));
    copyfile(fullfile(fn.fold_an_sess,[resBaseName '.eps']),rasters_folder);
    cellsArray(ic).resp=resp;
end

%save(fullfile(rasters_folder,'litralsArray.m'),'cellsArray');
load(fullfile(rasters_folder,'litralsArray.m'))
%get all latencies by concentration for every odor that has 3
%concentrations
observable = 'latency';
odors = {'2-hydroxyacetophenone','menthone'};
clear obs;
for io=1:numel(odors)
    odor=odors{io};
    concVector = [];
    obsMatrix  = nan(numel(cellsArray),3);
    for ic=1:numel(cellsArray)
        %get all cells that have 3 concentrations of the odor
        resp = cellsArray(ic).resp;
        stim=[resp.stim];
        theOdorAppeareances = find(strcmpi(odor,{stim.odorName}));
        if numel(theOdorAppeareances)<3
            continue;
        end
        
        resp=resp(theOdorAppeareances);
        %get the vector of observable
        if isempty(concVector)
            concVector = [stim(theOdorAppeareances).odorConc];
        end
        obsMatrix(ic,:) = [resp.(observable)];
    end
    obsMatrix(obsMatrix==0)=nan;
    obs(io).obs = observable;
    obs(io).odor=odor;
    obs(io).concVector = concVector;
    obs(io).obsMatrix = obsMatrix;
    theCell = cellsArray(ic);
    obs(io).rasterId = sprintf('%s_%03d_%s_odor_units%02d_resp',theCell.mouse,theCell.sess,theCell.rec,theCell.clu);
    obs(io).cellId = theCell.Id;
end

save(fullfile(rasters_folder,['litralsResp_' observable '.m']),'obs');
%plot me something nice
%first: just the latencies for each series:
%%
colors={'r','b'};
close all

io=1;
cv=obs(io).concVector;
cm=obs(io).obsMatrix;
h(io)=figure();
set(h(io),'Nextplot','add');
sh(io)=subplot(1,1,1);
suptitle(obs(io).odor);
hold on
for i=1:numel(cm(:,1))
    plot(sh(io),cm(i,:),i*[1 1 1],'-^',...
        'MarkerSize',16,'LineWidth',1,'Color',colors{io},...
        'MarkerFaceColor','w')
    text(cm(i,:),i*[1 1 1],{'1' '2' '3'},'VerticalAlignment','middle', ...
        'HorizontalAlignment','center')
end
ylabel(sh(io),sprintf('%s','cell'));
xlabel(sh(io),obs(io).obs);
set(sh(io),'XLim',[0 400])
set(sh(io),'YLim',[0 numel(cm(:,1))+1])
figName = sprintf('%s_conc_%s',obs(io).obs,obs(io).odor);
print(h(io), '-dpdf', fullfile(rasters_folder,[figName '.pdf']))





