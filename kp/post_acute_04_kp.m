function za=post_acute_04_kp(mouse,sess,rec)

za.get_info      = @get_info;
za.new_trials    = @za.new_trials;
za.file_names    = @file_names;
za.group_stimuli = @group_stimuli;
za.make_rasters  = @make_rasters;
za.first_analysis= @first_analysis;
za.structfind = @structfind;

end

function fn = file_names(mouse, sess, rec, stat)

%     local_disk  = '/home/zeke/data';
    local_disk='/Volumes/spikefolder';
%   local_disk  = '/experiment';
    stat_disk='/stations';
    %     server_disk = '';
    %     comp_disk   = '';

    fn.fold_rd_data  = fullfile(local_disk, 'raw_data');
    fn.fold_ss_data  = fullfile(local_disk, 'ss_data');
    fn.fold_pr_data  = fullfile(local_disk, 'pr_data');
    fn.fold_config   = fullfile(local_disk, 'SpikeGL_config');
    fn.fold_sd_data  = fullfile(stat_disk, '');


    if nargin > 0
        if isnumeric(mouse)
           mouse = sprintf('%04d', mouse);
        end
        if isnumeric(sess)
            sess = sprintf('%03d', sess);
        end

        fn.fold_rd_mouse = fullfile(fn.fold_rd_data,  sprintf('mouse_%s', mouse));
        fn.fold_rd_sess  = fullfile(fn.fold_rd_mouse, sprintf('sess_%s', sess));
        fn.log           = fullfile(fn.fold_rd_sess,  sprintf('log_%s_%s.txt', mouse, sess));

        fn.fold_ss_sess  = fullfile(fn.fold_ss_data,  sprintf('ss_%s_%s', mouse, sess));
        fn.fold_pr_sess  = fullfile(fn.fold_pr_data,  sprintf('m%s_%s', mouse, sess));
        fn.ss_sess_info  = fullfile(fn.fold_ss_sess,  sprintf('%s_%s_sess_info.mat', mouse, sess));
        fn.sess_info     = fullfile(fn.fold_pr_sess,  sprintf('%s_%s_sess_info.mat', mouse, sess));
    end

    if nargin >2 && ~isempty(rec)
        if isnumeric(rec)
            rec = sprintf('%02d', rec);
        end
        
        fn.fold_ss_rec  = fullfile(fn.fold_ss_sess,   sprintf('rec_%s',rec));
        fn.ss_rec       = fullfile(fn.fold_ss_rec,    sprintf('rec_%s.dat',rec));
        fn.rsm_data     = fullfile(fn.fold_pr_sess,   sprintf('%s_%s_%s_rsm.mat', mouse, sess, rec));
        fn.stim_type    = fullfile(fn.fold_pr_sess,   sprintf('%s_%s_%s_stim.mat', mouse, sess, rec));
        fn.raster       = fullfile(fn.fold_pr_sess,   sprintf('%s_%s_%s_raster.mat', mouse, sess, rec));
    end
    
    
    if nargin >3
        if isnumeric(stat)
            stat = sprintf('%02d', stat);
        end
        
        fn.fold_sd_data  = fullfile(stat_disk, sprintf('stat_%s',stat),'raw_data');
        fn.fold_sd_mouse = fullfile(fn.fold_sd_data,   sprintf('mouse_%s',mouse));
        fn.fold_sd_sess  = fullfile(fn.fold_sd_mouse,  sprintf('sess_%s',sess));
    end
end

function info=get_info(mouse,sess)
fn = file_names(mouse, sess)
load(fn.sess_info)
end

function stim=group_stimuli(mouse,sess,rec)
%gets the stimuli trhough all the runs
if ~exist('mouse', 'var')
    mouse='SDbehavingM72';
    sess='003';
    rec='e';
end

%get the files
fn=file_names(mouse, sess, rec);
%load trial
load(fullfile(fn.fold_pr_sess,sprintf('%s_%s_%s_trial.mat',mouse,sess,rec)));
%load unit;
load(fullfile(fn.fold_pr_sess,sprintf('%s_%s_%s_spikes.mat',mouse,sess,rec)));
%load sess_info
load(fn.sess_info)
rounding=2;

fprintf('*** Separating trial types for session %s, rec %s\n',sess,rec);
v=unique([trial.laserAmp]);
iStim=0;
for iV=1:numel(v)
    fprintf('V= %2.2f\n',v(iV))
    vTrials=find([trial.laserAmp]==v(iV));
    durs=unique(round(diff([trial(vTrials).laserTimes]*rounding))/rounding);
    for iDur=1:numel(durs)
        fprintf('  Dur= %2.1f\n',durs(iDur))
        iStim=iStim+1;
        stim(iStim).trials=find(round(diff([trial.laserTimes]*rounding))/rounding==durs(iDur) & [trial.laserAmp]==v(iV));
        stim(iStim).iD=iStim;
        stim(iStim).laserAmp=v(iV);
        stim(iStim).laserDur=durs(iDur);
        stim(iStim).laserPower=mean([trial(stim(iStim).trials).laserPower]);
        stim(iStim).runs=unique([trial(stim(iStim).trials).run]);
    end
end
fprintf('*** Found %d trial types for session %s, rec %s\n',iStim,sess,rec);
save(fn.stim_type,'stim');
end


%%
function [dataPoint] = make_rasters(mouse,sess,rec,cells,dp_U,lightral_cells,dp_L)
%gets the stimuli trhough all the runs
if ~exist('mouse', 'var')
    mouse='SDbehavingM72';
    sess='003';
    rec='e';
end

%get the files
fn=file_names(mouse, sess, rec);
%load trial
load(fullfile(fn.fold_pr_sess,sprintf('%s_%s_%s_trial.mat',mouse,sess,rec)))
%load unit;
load(fullfile(fn.fold_pr_sess,sprintf('%s_%s_%s_spikes.mat',mouse,sess,rec)))
%load sess_info
load(fn.sess_info)
%make or get the stimuli structure   
if ~exist(fn.stim_type, 'file')
    stim=group_stimuli(mouse,sess,rec);
else
    load(fn.stim_type)
end
% stim=group_stimuli(mouse,sess,rec);

%get the number of rec
if isnumeric(rec)
    rec = sprintf('%02d', rec);
end
iRec=structfind(info.rec,'name',rec)
Fs=info.rec(iRec).sampling_freq;


%% define units to inspect
if nargin < 4
    %for now only get good cells (lightral and mitral)
%     units = find([unit.qltMask]==7 | [unit.qltMask]==3)
else
    units = cells
    % units=[5 6]
    %for modification of clusters
    %lightral cells for ZKanesthTAAR4_002_a
    %clu(1).units=[3 7 8 9]
    if nargin < 5
%         lightral=find([unit.lightral])
        % lightral=find([unit.qltMask]==7)
    else
        lightral = lightral_cells;
    end
end

clu(1).units=lightral;
clu(1).dp=dp_L;   % save a datapoint matrix or not

%reomve units already in clusters;
for iClu=1:numel(clu)
    for iUnitClu =1:numel(clu(iClu).units)
        units( find( units==clu(iClu).units(iUnitClu) ) )  =[]
    end
end

%the rest are clusters as they come
for iUnits=1:numel(units)
    iClu=iClu+1;
    clu(iClu).units=units(1)
    clu(iClu).dp=dp_U;
    units(1)=[];
end
%%

%for every unit:
%get the rasters for all the stimuli
%plot them
bSpan=500;
aSpan=1000;
psthBin=1;
smoothWindow=2;

%segment to fit the pdf to the psth
t0Pdf=0;
tfPdf=30;


%% 
%for all the units in each cluster, get the raster, and fit a normal
%distribution to the segment between 0 and tfPdf (in miliseconds)

for iClu=1:numel(clu);
    thisClu=clu(iClu)
    clear figure(iClu)
    dpCount=0;
    
    for iStim=1:numel(stim)
        raster=[];
        ir=0;
        trialSelect=stim(iStim).trials;
        unitsReport=[];
        %psth originally in millisecods
        psth=zeros(1,round((aSpan+bSpan)/psthBin));
        psth_aux=zeros(size(psth));
        tPsth=(1:ceil(length(psth)/smoothWindow))*psthBin*smoothWindow-round(bSpan)-smoothWindow/2;
        %gather the rasters for all the units with this stimulus in this
        %cluster
        for iUnit=1:numel(thisClu.units)
            unitId=thisClu.units(iUnit);
            unitsReport=[unitsReport num2str(unitId) '-'];
            for it=1:numel(trialSelect)
                trN=trialSelect(it);
                iRun=trial(trN).run;
                runOffset=info.rec(iRec).run(iRun).offset/Fs*1000;
                t0=trial(trN).laserTimes(1)+trial(trN).start+runOffset;
                tspikes=unit(unitId).times;
                if numel(raster)<it
                    raster(it).t=[];
                end
                raster(it).t=[raster(it).t ;tspikes(find(tspikes > (t0 - bSpan) & tspikes < (t0 + aSpan))) - t0];
                %set the segment of stimulation to zero (remove those
                %spikes that are likely artifacts)
                raster(it).t( 0<raster(it).t & raster(it).t<diff(trial(trN).laserTimes)+0.1)=[];
                psth_aux(:)=0;
                psth_aux(ceil(raster(it).t/psthBin+bSpan/psthBin))=1;
                psth=psth+psth_aux;
            end
        end
        figure(iClu)
        gs=subplot(ceil(numel(stim)/3),3,iStim);
        hold on
        for ir=1:length(raster)
            plot(raster(ir).t,ones(size(raster(ir).t))+ir*1,'k.','MarkerSize',5);
        end
        title(sprintf('%2.2f V (%2.2f mW), %2.1f ms',stim(iStim).laserAmp,stim(iStim).laserPower,stim(iStim).laserDur))
        set(gs,'FontSize',6,'Xlim',[-50,100]);
        box off
        psthSmooth=smooth(psth,smoothWindow);
        psthSmooth=resample(psth,1,smoothWindow);
        %remove the baseline and normalize
        psthSmooth=(psthSmooth-mean(psthSmooth(1:round(bSpan/(psthBin*smoothWindow)))))/numel(psthSmooth);
        plot(tPsth,psthSmooth*numel(raster),'b');
        
        
        %% 
        %make a poisson distribution out of the post-stimulus part of the psth
        
        % t0Poi=diff(trial(trN).laserTime);
        pdfSelect=tPsth>t0Pdf & tPsth<tfPdf;
        psthPdf=psthSmooth(pdfSelect);
%         totalSpikes=mean(psth(t0Pdf*psthBin +bSpan : tfPdf*psthBin+bSpan))-mean(psth(tPsth<0))
        baseline_window = tPsth<0;
        baseline_spikes_sec = sum(psth(baseline_window))/length(baseline_window)*1000/psthBin/smoothWindow/numel(trialSelect);
        resp_window = t0Pdf*psthBin +bSpan : tfPdf*psthBin+bSpan;
        resp_spikes_sec = sum(psth(resp_window))/length(resp_window)*1000/psthBin/smoothWindow/numel(trialSelect);
        avgChange_spikes_sec = resp_spikes_sec - baseline_spikes_sec;
        psthPdf(psthPdf<0)=0;
        tPdf=tPsth(pdfSelect);
        %normalize the psth in the range chosen for distribution
        psthPdf=psthPdf/sum(psthPdf.*tPdf);
        plot(tPdf,psthPdf/max(psthPdf)*numel(trialSelect),'r')
        
        %if it has a significant peak, fit a distribution
        if(max(psthPdf)>2.5*std(psthPdf))
            %get all the times of times of all the spikes
            tspike=[];
            firstSpike=[];
            for ir=1:numel(raster)
                rasterCut=transpose( raster(ir).t(raster(ir).t>t0Pdf & raster(ir).t<tfPdf) );
                if ~isempty(rasterCut)
                    tspike=[tspike rasterCut];
                    firstSpike=[firstSpike rasterCut(1)];
                else
                    firstSpike=[firstSpike 0];
                end
            end
            [mu, sigma]=normfit(tspike)
            if (sigma)
                spikeDist=pdf('normal',tPdf,mu,sigma)/sum(tPdf.*pdf('normal',tPdf,mu,sigma));
                plot(tPdf,(spikeDist)/max(spikeDist)*numel(trialSelect),'g')
            end
        end
      %%
        
        %save the responses in a structure
        %stimulus, cluster, response
        if clu(iClu).dp
            dpCount=dpCount+1;
            dpAux.stim=stim(iStim);
            dpAux.stimV=stim(iStim).laserAmp;
            dpAux.stimDur=stim(iStim).laserDur;
            dpAux.stimPow=stim(iStim).laserPower;
            dpAux.stimEner=stim(iStim).laserAmp*stim(iStim).laserDur;
            dpAux.units=thisClu.units;
            dpAux.raster=raster;
            dpAux.psth=[tPsth;psthSmooth];
            dpAux.avgSpikes=avgChange_spikes_sec;
            dpAux.delay=mean(firstSpike(firstSpike>0));
            dpAux.reliable=sum(firstSpike>0)/numel(firstSpike);  % proportion of trials that had a spk in 50ms from laser off
            if(max(psthPdf)>2.5*std(psthPdf))
                dpAux.pdf='normal';
                dpAux.pdfMu=mu;
                dpAux.pdfSigma=sigma;
            else
                dpAux.pdf= '';
                dpAux.pdfMu=nan;
                dpAux.pdfSigma=nan;
            end
            dataPoint(dpCount)=dpAux;
            clear dpAux
            save(fullfile(fn.fold_pr_sess,sprintf('%s_%s_%s_%d_datapoints.mat',mouse,sess,rec,3)),'dataPoint');
        end
    end
    unitsReport(end)=[];
%     disp(unitsReport)
    suptitle(['Mouse ' mouse ', rec ' sess rec '. Unit ' unitsReport]);
    print('-dpdf',fullfile(fn.fold_pr_sess,sprintf('%s_%s_%s_%s_raster.pdf',mouse,sess,rec,unitsReport)));
end

fprintf('*** Found %d trial types for session %s, rec %s\n',iStim,sess,rec);
end


%%

function [vf,df,dataPoint]=first_analysis(mouse,sess,rec,iClu)
%assume cluster (of cells in datapoints structure)=1;
if ~exist('mouse', 'var')
    mouse='SDbehavingM72';
    sess='003';
    rec='e';
    iClu=1;
end

%get the files
fn=file_names(mouse, sess, rec);
%load datapoints structure
dp.name=fullfile(fn.fold_pr_sess,sprintf('%s_%s_%s_%d_datapoints.mat',mouse,sess,rec,iClu));
load(dp.name);

vstims=unique([dataPoint.stimV])
dstims=unique([dataPoint.stimDur])

%plot a property grouped by V
%prop='avgSpikes';
props={'pdfSigma','pdfMu','delay','avgSpikes'} ;
vf=figure(1)
hold on;
for iprop=1:numel(props)
prop=props{iprop};
toplot=zeros(length(dstims),length(vstims))
for iv=1:length(vstims)
    v=vstims(iv)
    leg{iv}=[num2str(v) 'V'];
    subset=dataPoint([dataPoint.stimV]==v)
    dselect=[subset.stimDur]
    full=find(ismember(dstims,dselect))
    empty=find(~ismember(dstims,dselect))
    for ifull=1:numel(full)
        toplot(full(ifull),iv)=subset(ifull).(prop)
    end
    for iempty=1:numel(empty)
        toplot(empty(iempty),iv)=nan;
    end
end
subplot(2,2,iprop)
plot(dstims,toplot,'-*')
xlabel('duration (ms)');
ylabel(sprintf('%s',prop));
legend(leg)
end
suptitle(['Mouse ' mouse ', rec ' sess rec '. Units group ' num2str(iClu)]);
fnFig=fullfile(fn.fold_pr_sess,sprintf('%s_%s_%s_%d_d.pdf',mouse,sess,rec,iClu));
print('-dpdf',fnFig);
hold off

%plot a property grouped by duration
df=figure(2)
hold on
for iprop=1:numel(props)
prop=props{iprop};
toplot=transpose(zeros(length(dstims),length(vstims)))
for id=1:length(dstims)
    d=dstims(id)
    leg2{id}=[num2str(d) 'ms'];
    subset=dataPoint([dataPoint.stimDur]==d)
    vselect=[subset.stimV]
    full=find(ismember(vstims,vselect))
    empty=find(~ismember(vstims,vselect))
    for ifull=1:numel(full)
        toplot(full(ifull),id)=subset(ifull).(prop)
    end
    for iempty=1:numel(empty)
        toplot(empty(iempty),id)=nan;
    end
end
subplot(2,2,iprop)
tf=plot(vstims,toplot,'-^')
tleg=legend(leg2);
xlabel('V stim (V)');
ylabel(sprintf('%s',prop));

end
suptitle(['Mouse ' mouse ', rec ' sess rec '. Units group ' num2str(iClu)]);
fnFig=fullfile(fn.fold_pr_sess,sprintf('%s_%s_%s_%d_v.pdf',mouse,sess,rec,iClu));
print('-dpdf',fnFig);
hold off

figure; 
scatter([dataPoint.StimEner],[dataPoint.reliable])
xlabel('Laser stim energy')
ylabel('Reliability')

end
function index=structfind(a,field,value)
% StructFind, Find the index of a certain string or value in a struct
%
%       index=structfind(a,field,value)
%
%  inputs,
%       a : A Matlab struct, for example a(1).name='red', a(2).name='blue';
%       field : The name of the field which is searched, for example 'name'
%       value : The search value, for example 'blue'
%
%  outputs,
%       index : The Struct index which match the search
%
%

% We don't compare structs
if(isstruct(value)),
    error('structfind:inputs','search value can not be a struct');
end

% Stop if field doesn't exist
if(~isfield(a,field))
    index=find(arrayfun(@(x)(cmp(x,field,value)),a,'uniformoutput',true));
else
    index=find(arrayfun(@(x)(cmp(x,field,value)),a,'uniformoutput',true));
end


    function check=cmp(x,field,value)
        check=false;
        if(isfield(x,field))
            % Simple field like x.tag
            x=x.(field);
        else
            % Complex field like x.tag.child.value
            in=find(field=='.');
            s=[1 in+1]; e=[in-1 length(field)];
            for i=1:length(s)
                fieldt=field(s(i):e(i));
                if(isfield(x,fieldt)), x=x.(fieldt);  else return; end
            end
        end
        
        % We don't compare structs
        if(isstruct(x)), return; end
        
        % Values can only be equal, if they equal in length
        if(length(x)==length(value)),
            % This part compares the NaN values
            if((~iscell(x))&&(~iscell(value))&&any(isnan(value))),
                checkv=isnan(value); checkx=isnan(x);
                if(~all(checkx==checkv)), return; end
                x(checkx)=0; value(checkv)=0;
            end
            % This part compares for both string as numerical values
            if(iscell(x)||iscell(value))
                check=all(strcmp(x,value));
            else
                check=all(x==value);
            end
        end
    end
end


function index=structfind(a,field,value)
% StructFind, Find the index of a certain string or value in a struct
%
%       index=structfind(a,field,value)
%
%  inputs,
%       a : A Matlab struct, for example a(1).name='red', a(2).name='blue';
%       field : The name of the field which is searched, for example 'name'
%       value : The search value, for example 'blue'
%
%  outputs,
%       index : The Struct index which match the search
%
%

% We don't compare structs
if(isstruct(value)),
    error('structfind:inputs','search value can not be a struct');
end

% Stop if field doesn't exist
if(~isfield(a,field))
    index=find(arrayfun(@(x)(cmp(x,field,value)),a,'uniformoutput',true));
else
    index=find(arrayfun(@(x)(cmp(x,field,value)),a,'uniformoutput',true));
end


    function check=cmp(x,field,value)
        check=false;
        if(isfield(x,field))
            % Simple field like x.tag
            x=x.(field);
        else
            % Complex field like x.tag.child.value
            in=find(field=='.');
            s=[1 in+1]; e=[in-1 length(field)];
            for i=1:length(s)
                fieldt=field(s(i):e(i));
                if(isfield(x,fieldt)), x=x.(fieldt);  else return; end
            end
        end
        
        % We don't compare structs
        if(isstruct(x)), return; end
        
        % Values can only be equal, if they equal in length
        if(length(x)==length(value)),
            % This part compares the NaN values
            if((~iscell(x))&&(~iscell(value))&&any(isnan(value))),
                checkv=isnan(value); checkx=isnan(x);
                if(~all(checkx==checkv)), return; end
                x(checkx)=0; value(checkv)=0;
            end
            % This part compares for both string as numerical values
            if(iscell(x)||iscell(value))
                check=all(strcmp(x,value));
            else
                check=all(x==value);
            end
        end
    end
end