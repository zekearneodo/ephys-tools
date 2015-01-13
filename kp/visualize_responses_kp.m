% Script for visualizing spikes aligned with the onset of stimuli
% uses data_managment_tools_03
% inits the parameters and the functions for analysis,
% leaves global vp with the parameters for all the functions to call
%has functions:
% -find_stimuli_set

function vr = visualize_responses_kp(mouse,sess,rec,unit,sType)
global vp vd dm;  

%add the package of common files and folders to the path.
%here are the functions for file_naming, for instance, and any other
%functions that several script bundles will use.
%that is /baseFolder/ephysDataManagement/current/include
includePath=fullfile(fileparts(pwd),'current','include');
addpath(includePath);
% Functions of include that it uses:
%   - file_names(mouse,sess,rec,stat)

vr.visual_par_init     = @visual_par_init;
vr.make_stimuli_set    = @make_stimuli_set;
vr.view_stimuli_set    = @view_stimuli_set;
vr.analysis_file_names = @analysis_file_names;
% pp.ss_prep        = @ss_prep;
% pp.resampling     = @resampling;

if nargin < 5
    sType='laser'
end

dm = data_management_tools_032();
vp = visual_par_init(mouse,sess,rec,unit,sType);
vd = make_stimuli_set(vp);

view_stimuli_set(sType,unit);

end

function vp=visual_par_init(mouse,sess,rec,unit,sType)
% visual_par_init(mouse,sess,rec,[unit])
%initializes or resets the parameters for the visualization
%that will be stored in the global vp;
%makes the stimuli set
global dm vd;

vp.mouse = mouse;
vp.sess  = sess;
vp.rec   = rec;

%un    = 13;

%     fn.fold_sd_data  = fullfile(stat_disk, '');
if nargin > 1
    if isnumeric(mouse)
        vp.mouse = sprintf('%04d', mouse);
    end
    if isnumeric(sess)
        vp.sess = sprintf('%03d', sess);
    end
end

if nargin >2 && ~isempty(rec)
    if isnumeric(rec)
        vp.rec = sprintf('%02d', rec);
    end
end

if strcmp(sType,'laser')
    vp.t1=-200;
    vp.t2=400;
else
    vp.t1  = -200;
    vp.t2  = 400;
end
    
vp.nt  = vp.t2-vp.t1;
vp.bin = 15;
vp.t   = mean(reshape((vp.t1+1):vp.t2, vp.bin , vp.nt/vp.bin),1);

vp.par = {'odorName', 'odorConc', 'laserDur', 'laserAmp', 'pulseGroup','pulseOffset'};
vp.stimTypes={'laser','odor'};

vp.fn = dm.file_names(mouse, sess, rec);
vp.fn = analysis_file_names(mouse,sess,rec);

q = load(vp.fn.trials);

%check if there is a trial correction function to be applied
trCFname=fullfile(vp.fn.fold_pr_sess,sprintf('%s_%s_trial_correct.m', vp.mouse,vp.sess));
if exist(trCFname)
    fprintf('Found trial correction function %s\n',trCFname);
    p = path;
    path(p,vp.fn.fold_pr_sess);
    fun=sprintf('%s_%s_trial_correct', vp.mouse,vp.sess);
    trCorrectFcn=eval(['@' fun ';'])
    path(p);
    trial=feval(trCorrectFcn,q.trials);
    vp.tr=trial;
else
vp.tr = q.trial;
end

%% select trials to look at

% % find indices of et tig trials
% for ii = 1:numel(vp.tr)
%     odorName{ii} = vp.tr(ii).odorName;
% end
% et = strfind(odorName,'ethyl tiglate');
% jj = 1;
% for ii = 1:numel(vp.tr)
%     if et{ii} == 1
%         et_in(jj) = ii;
%         jj = jj+1;
%     end
% end

% % find first trial after certain time in ms
% first_tr = find([vp.tr.start]>600000,1);

% select the trials you want to keep;
% in_tr=(first_tr:numel(vp.tr));
% vp.tr=vp.tr(in_tr)



vd=make_stimuli_set(vp);
end

function fn=analysis_file_names(mouse,sess,rec)
% Load the filenames of the data, as in data_managment_tools,
% and create folders to store the fist stages of analysis
%add the package of common files and folders to the path.
%here are the functions for file_naming, for instance, and any other
%functions that several script bundles will use.
%that is /baseFolder/ephysDataManagement/current/include
includePath=fullfile(fileparts(pwd),'current','include');
addpath(includePath);
% Functions of include that it uses:
%   - file_names(mouse,sess,rec,stat)

global dm;

if isnumeric(mouse)
    mouse = sprintf('%04d', mouse);
end
if isnumeric(sess)
    sess = sprintf('%03d', sess);
end
if isnumeric(rec)
    rec = sprintf('%02d', rec);
end
        
fn = dm.file_names(mouse, sess, rec);

fn.fold_an_mouse  = fullfile(fn.disk,'analysis',sprintf('mouse_%s',mouse));
fn.fold_an_sess   = fullfile(fn.fold_an_mouse,sprintf('sess_%s',sess));
fn.basename_an    = sprintf('%s_%s_%s_',mouse,sess,rec);
end

function vd=make_stimuli_set(vp)
% ====================================================================================
% finding stimulus set, and filtering by properties of stimulation
% eg: odorstim,lightstim

stim = struct('odorName', '', 'odorConc', 0, 'laserDur', 0, 'laserAmp', 0, 'in_tr', [], 'pulseGroup', [],'pulseOffset',[])

for it=1:numel(vp.stimTypes)
    sType=vp.stimTypes{it}
    vd.(sType) = struct('list',[],'sort','')
end


%% Group the laser pulse times according to sniff cycle
% min is 21
ms_from_inh = [ 21 28;   68 75;   77 78;   101 108];  %these are the windows in ms!

    for jt = 1:numel(vp.tr)
        vp.tr(jt).pulseGroup = NaN;
    end
    used = [];
    for grp = 1:size(ms_from_inh,1)
        grp_ind = find([vp.tr.pulseOffset]>ms_from_inh(grp,1) & [vp.tr.pulseOffset]<ms_from_inh(grp,2));
        used = [used grp_ind];
        for jt = 1:numel(vp.tr)
            if ismember(jt,grp_ind)
                vp.tr(jt).pulseGroup = grp;
                vp.tr(jt).pulseOffset=round(mean(ms_from_inh(grp,:)));
            end
        end
    end
%%
ns   = 1;
tr=vp.tr(find([vp.tr.pulseGroup]>0));
par=vp.par;

for it = 1:numel(tr)
    if tr(it).start == 0
        continue
    end
        
%     tr(it).laserDur = diff(tr(it).laserTimes);
    tr(it).laserDur =1;
    new = 1; ks = 0;
    while new&&(ks<ns)
        ks = ks + 1;
        same_stim = strcmp(tr(it).odorName, stim(ks).odorName);
        for ip = 2:numel(par)
            same_stim = same_stim&(tr(it).(par{ip}) == stim(ks).(par{ip}));
        end
        
        if same_stim
             new = 0;
             stim(ks).in_tr(end+1) = it;
        end
    end
    
    if new
        ns = ns + 1;
        for ip = 1:numel(par)
            stim(ns).(par{ip}) = tr(it).(par{ip});
        end
        stim(ns).in_tr = it;
    end
end

stim = stim(2:end)
ns   = ns - 1;


%remove from the list the stimuli for which there are less than 10 trials
keepStim=[];
for is=1:numel(stim)
    if numel(stim(is).in_tr)>10
        keepStim=[keepStim is];
    end
end
stim=stim(keepStim);


%select stimuli types
for it=1:numel(vp.stimTypes)
    sType=vp.stimTypes{it};
    switch sType
        case 'odor'
            sSelect = find([stim.odorConc]>0);
            sSort   = 'odorName';
        case 'laser'
            sSelect = find([stim.odorConc]==0);
%             sSort   = 'laserAmp';
            sSort   = 'pulseGroup';
    end
    vd.(sType).list = sSelect;
    vd.(sType).sort = sSort;
end
vd.stim=stim;
end %function vd=make_stimuli_set()


function view_stimuli_set(sType,un)
global vd vp;

stim=vd.stim(vd.(sType).list);
ns=numel(stim)

sSort=vd.(sType).sort;

if isnumeric(stim(1).(sSort))
    [~, in_sort] = sort(-[stim.(sSort)]);
else
    [~, in_sort] = sort({stim.(sSort)});
end

stim = stim(in_sort);
%selection of stimuli for visualization of mouse 1461, sess 004, reca
%stim14614a=[6 1 5 7]
%stim = stim(stim14614a);

%in_stim = [1,3,8,6,7];
% stim(1).in_tr = stim(1).in_tr(51:300);
% stim(3).in_tr = stim(3).in_tr(10:55);

% stim(1).lineColor = 'k';
% stim(2).lineColor = 'k';
% stim(3).lineColor = 'b';
% stim(4).lineColor = 'k';
ns=numel(stim);
for i=1:ns
    stim(i).lineColor = 'r';
end

for i=1:ns
    stim(i).lineStyle = '-';
    stim(i).lineColor = 'k';
    if stim(i).odorConc < 0.014
        stim(i).lineStyle = '-';
        stim(i).lineColor = 'r';
    end
end

% stim(1).lineStyle = '-';
% stim(3).lineStyle = '-';
% stim(8).lineStyle = '--';
% stim(6).lineStyle = '-';
% stim(7).lineStyle = '--';
rmax = 200;

% ===============================================
% plotting raster plots and psth for a group of units
ifigRas=12;
ifigPsth=22;
gf = figure(ifigRas); clf
ip = 0;
tr=vp.tr;
t1=vp.t1;
t2=vp.t2;
bin=vp.bin;
t=vp.t;
nt=vp.nt;
fn=vp.fn;

% quick get baseline for all the odor trials, two sniffs before final valve
% opening
if strcmp(sType,'odor')
    rateBase=zeros(1,vp.nt);
    nb=0;
    for kt=1:numel(tr)
        trial=tr(kt);
        if isempty(trial.sniffZeroTimes) || isempty(trial.spikeTimes)
            continue
        end
        iiBase = find(trial.sniffZeroTimes(1,:)>0,1) -3;
        if isempty(iiBase) | iiBase<1
            continue
        end
        sp = sort(round(vertcat(trial.spikeTimes{un})));
        t0Base = trial.sniffZeroTimes(1,iiBase)
        spBase = sp((sp>t0Base+t1)&(sp<t0Base+t2))-t0Base;
        rateBase(spBase-t1)= rateBase(spBase-t1)+1;
        nb=nb+1;
    end
    rateBase = sum(reshape(rateBase, bin, nt/bin),1)/kt/(bin/1000);
end
    

for ks = 1:numel(stim)
%     rmax=numel(stim(ks).in_tr);
    rate = zeros(1,vp.nt);
    x    = zeros(1e5,1);    y    = zeros(1e5,1);
    nsp  = 0;
    kt   = 0;

    upperTrials=min(20,numel(stim(ks).in_tr));
    stim(ks).in_tr=stim(ks).in_tr(end-upperTrials+1:end);
    spkCount=[];
    for it = stim(ks).in_tr
        %get the spikeTimes for all the units
        
        sp = sort(round(vertcat(tr(it).spikeTimes{un})));

        %if its odor align with sniff, otherwise align with laser
        %presentation
        if stim(ks).odorConc >0 && strcmp(sType,'odor')
            %find first sniff after final valve open
            ii = find(tr(it).sniffZeroTimes(1,:)>tr(it).odorTimes(1)*1.009,1)
            if isempty(ii) || ii<1
                continue
            end
            t0 = tr(it).sniffZeroTimes(1,ii)
            stim_str = sprintf('%s \n %1.1e', tr(it).odorName(1:min(length(tr(it).odorName),10)), tr(it).odorConc);
        
        else
            t0 = tr(it).laserTimes(1);
            stim_str = sprintf('Laser: \n%2.1fV -%3.1fms; Inh delay: %i', stim(ks).laserAmp/1000, stim(ks).laserDur, stim(ks).pulseGroup);
        end
        kt = kt + 1;
        
        sp = sp((sp>t0+t1)&(sp<t0+t2))-t0;
        rate(sp-t1) = rate(sp-t1) + 1;
        n  = numel(sp);

        %another way of counting the number of spikes every trial
        %(how many spikes are from t0 to t0+200); but t0 is now the 0 in sp
        spkCount=[spkCount numel(sp((sp>0)&(sp<200)))];
        x(nsp+(1:n)) = sp(:);   
        y(nsp+(1:n)) = kt*ones(1,n);
        nsp = nsp + n;
        
   
    end
    x = x(1:nsp);
    y = y(1:nsp);
    rate = sum(reshape(rate, bin, nt/bin),1)/kt/(bin/1000);
%     rateBase = sum(reshape(rateBase, bin, nt/bin),1)/kt*100;
       
   
    baseline=mean(rate(t<0));
    pkThresh=baseline + 3*std(rate(t<0));
    [pkRate,pkInd]=findpeaks(rate(t>0),'MINPEAKHEIGHT',pkThresh,'NPEAKS',1);
    pkTime=t(pkInd+find(t>0,1)-1);
    spkBase=sum(rate(t<0))*(bin/1000)
    
    if isempty(pkTime)
        pkTime=0;
        pkRate=0;
    end
    nSpikes=(sum(rate(t>0 & t<200)*bin))/1000
    vd.resp(ks).spkCount  = spkCount;
    vd.resp(ks).spkBase   = spkBase;
    vd.resp(ks).ntrial    = numel(stim(ks).in_tr);
    vd.resp(ks).nSpikes   = nSpikes;
    vd.resp(ks).nSpk      = mean(spkCount);
    vd.resp(ks).nSpkErr   = std(spkCount);
    vd.resp(ks).delay     = pkTime;
    vd.resp(ks).maxFR     = pkRate;
    vd.resp(ks).baseline  = baseline;
    vd.resp(ks).stim      = stim(ks);

    %%%
 
     %plots this psth
    ip = ip + 1;
    gf = figure(ifigRas)
    gs = subplot(2,numel(stim), ip);
    yt  = max([1, floor(kt/10)*10]);
    plot(x,y, '.', 'MarkerSize',7), hold on
    plot([0,0], [0, kt+1], '--k'), hold off
    set(gs, 'XLim', [t1, t2], ...
        'YLim',      [0, kt+1], ...
        'YTick',     [0, yt], ...
        'XTick',     [-200,0,200], ...
        'FontSize',  10);
    title(stim_str, 'FontSize', 10, 'FontWeight', 'bold')

    
    gs = subplot(2,numel(stim), ip+numel(stim));
    plot(t, rate, ...
        'LineWidth', 1, ...
        'Color',     stim(ks).lineColor, ...
        'LineStyle', stim(ks).lineStyle),
    hold on
    %plot(t,ones(size(t))*baseline,'.k')
    %plot(t,ones(size(t))*pkThresh,'b')
    if strcmp(sType,'odor')
        plot(t,rateBase,'LineStyle','--','Color',[.5,0.5,0.5])
    end
    plot(pkTime,pkRate,'k*')
    plot([0,0], [0, rmax], '--k'), hold off
    set(gs, 'XLim',     [t1, t2], ...
        'XTick',         [-200,0,200], ...
        'YLim',          [0, rmax], ...
        'YTick',         0:20:rmax, ...
        'FontSize',      10);
    
    %plots of psths altogether
    if ks>1
        figure(ifigPsth)
        gs0 = subplot(1,1,1);
        plot(t, rate, ...
            'LineWidth', 2, ...
            'Color',     stim(ks).lineColor, ...
            'LineStyle', stim(ks).lineStyle),
        hold on
    end
    
    
end

gf = figure(ifigRas)
sUnits=num2str(un(1));
for iu=2:numel(un)
    sUnits=sprintf('%s,%2d',sUnits,un(iu));
end
suptitle(sprintf('Mouse %s - session %3d - rec %s units %s',num2str(vp.mouse),str2num(vp.sess),vp.rec,sUnits));
%saves the response data (the plot and the response structure)
if ~exist(fn.fold_an_mouse, 'dir')
    mkdir(fn.fold_an_mouse)
end
if ~exist(fn.fold_an_sess, 'dir')
    mkdir(fn.fold_an_sess)
end

fileBase=[fn.basename_an sType '_units'];
for iu=1:numel(un)
    fileBase=[fileBase num2str(un(iu),'%02d')]
end

save(fullfile(fn.fold_an_sess,sprintf('%s_resp.mat',fileBase)), '-struct', 'vd','resp')
print(gf,'-depsc',fullfile(fullfile(fn.fold_an_sess,sprintf('%s_resp.eps',fileBase))))
    

if ks>1
    figure(ifigPsth)
    gs0 = subplot(1,1,1);
    subplot(gs0)
    plot([0,0], [0, rmax], '--k'), hold off
    set(gs0, 'XLim',     [t1, t2], ...
        'XTick',         -200:100:200, ...
        'YLim',          [0, rmax], ...
        'YTick',         0:20:rmax, ...
        'FontSize',      10);
end
end
