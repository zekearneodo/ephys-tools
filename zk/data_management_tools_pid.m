%PID data_managment tools
%01-26-2014
%CW (core) / ZK (adaptation for dilluter and modifications)

%This set of functions reads h5 files generated by Voyeur with the pid
%trace recorded in the sniff channel, and performs some anaylyses.
%Assumes the file structure is the same as for mice (pid are taken as
%mice).
%The only difference is that it look for x_yy_*_pid.h5 files (x rec, yy
%run).

% -pid_analysis(mouse,sess,stat) gets the pid traces for all the trials and does an analysis
%of the traces and the mean values as a function of odor concentration for
%all the configurations (vialConc, dillution) of all odors.
%It pulls the data if stat is specified, otherwise assumes it's already on
%the server.

% -pid_all_series(trPid) does the same analysis as above but it takes a
% trial structure with pre-selected conditions (odor names, for instance)

%-Pid_get(mouse,sess,rec) gets the h5 files for a pid, sess, all the recs,
%and leaves it into a trial structure trPid.

%-[concPlot]=pid_multi_graph(tr_pid)(pidStruct,odor,*odorConc,*vialConc,*dillution,suppress_graphs)
% plots traces and mean values of pid activity for odor (string), as a
% function of the vector odorConc[], grouped by conditions of vialConc[]
% and dillution[]
% retunrs a matrix with amplitudes, vector of concentrations for the x
% axis, and vectors of legends (what's in every row of amplitude matrix).


function dm = data_management_tools_pid()
global dm

% include the folder where the common files are usually located
[~,computerName]=system('hostname');
if ~strcmp(strtrim(computerName),'flipper')
    includePath=fullfile(fileparts(pwd),'current','include');
    addpath(includePath);
end

    dm.pid_get            = @pid_get;
    dm.pull_data          = @pull_data;
    dm.pid_graph          = @pid_graph;
    dm.pid_multi_graph    = @pid_multi_graph;
    dm.pid_max_graph      = @pid_max_graph;
    dm.pid_graph_all      = @pid_graph_all;
    dm.pid_all_series     = @pid_all_series;
    dm.pid_analysis       = @pid_analysis;
    dm.subplot_fig        = @subplot_fig;
    dm.vial_lookup        = @vial_lookup;
    
end

function pull_data(mouse, sess, stat)
% gets the data from mouse and session from the recording station stat
    fprintf('Pulling data from station %s \n',stat);
    fprintf('============================================================\n')
    fn  = file_names(mouse, sess,'',stat)
    
    %check existence of folder in the station, and if it is mounted
    if ~exist(fn.fold_sd_data, 'dir')
        error(['no raw data folder (raw_data) in station ' stat '. Is it mounted?'])
    end
    
    if ~exist(fn.fold_sd_mouse, 'dir')
        error(['no mouse ' mouse ' in station ' stat])
    end
    
    if ~exist(fn.fold_sd_sess, 'dir')
        error(['no session ' sess ' in station ' stat]);
    end
    
    %check target folder (raw data)

    if ~exist(fn.fold_rd_mouse, 'dir')
        mkdir(fn.fold_rd_mouse)
    end
    
    if ~exist(fn.fold_rd_sess, 'dir')
        mkdir(fn.fold_rd_sess)
    end
    
    % copy files one by one
    fl=dir(fn.fold_sd_sess);
    if numel(fl)<3
        error(['No files in ' fn.fold_sd_sess]);
    end
    
    isub=[fl(:).isdir];
    nameFolds = {fl(isub).name}';
    fl(ismember(nameFolds,{'.','..'})) = [];
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
    
    fprintf('Copying %d files (incl %d folders) from %s  \n to %s \n',numel(fl),numel(nameFolds),fn.fold_sd_sess,fn.fold_rd_sess);
    fprintf('===========================================================================\n')
    tic;
    
    for ifl = 1:numel(fl)
        source=fullfile(fn.fold_sd_sess,fl(ifl).name);
        dest=fn.fold_rd_sess;
        FileOrFolder='File';
        if(fl(ifl).isdir)
            dest=fullfile(fn.fold_rd_sess,fl(ifl).name);
            if ~exist(dest,'dir')
                mkdir(dest);
            end
            FileOrFolder='*(Folder)';
        end
        fprintf('%s %2d:\t %24s %5d Mb...',FileOrFolder,ifl,fl(ifl).name,round(fl(ifl).bytes/1000000));
        [status]=copyfile(source,dest);
        if ~status
           error(['Error copying ' source]);
        end
        fprintf(' ok \t %5d \n',round(toc));
    end
    
    
    
    fprintf('Done fetching files. \n');
    fprintf('===========================================================================\n')
end

function trPid=pid_get(mouse,sess,rec)
% This function reads an h5 file with PID traces recorded in the sniff
% channel.
%it returns a trial structure with the fields defined in var:
%new: is the name that it will have in the trPid struct.
%old: is the name to look for for the field in the h5, voyeur generated
%file.
%There are two fields that need to be read from the h5 file: the trial
%start and the trial end. the variables nTrEndField and nTrStart field
%refer to the index in the var() structure where those fields are defined
%(the program will look for their 'old' value when it needs to read them in
%the table.
%the new names can be whichever, and so do the old ones (names in
%global tr;
    if isnumeric(mouse)
        mouse=sprintf('%04d',mouse)
    end

    var(1)  = struct('new', 'trialNum',   'old', 'trialNumber', 'relative', 0);
    var(2)  = struct('new', 'trialstart', 'old', 'starttrial',  'relative', 0);
    var(3)  = struct('new', 'trialEnd',   'old', 'endtrial',    'relative', 0);
    var(4)  = struct('new', 'fvOnTime',   'old', 'fvOnTime',    'relative', 1);
    var(5)  = struct('new', 'fvdur',      'old', 'fvdur',       'relative', 0);
    
    var(6)  = struct('new', 'odor',       'old', 'odor',        'relative', 0);
    var(7)  = struct('new', 'vial',       'old', 'vial',        'relative', 0);
    var(8)  = struct('new', 'vialConc',   'old', 'vialconc',    'relative', 0);
    var(9)  = struct('new', 'air',        'old', 'AirFlow_1',   'relative', 0);
    var(10) = struct('new', 'nitrogen',   'old', 'NitrogenFlow_1','relative', 0);   
    var(11) = struct('new', 'odorConc',   'old', 'odorconc',    'relative', 0);
    var(12) = struct('new', 'dillution',  'old', 'dillution',   'relative', 0);
    
    % number of fields referring to trial start and end in the variables
    % description var() struct:
    nTrEndField = 3;
    nTrStartField=2;
    
    brEvents  = { 'breaks'};
    streams = {'sniff'};
   
    if isnumeric(mouse)
        mouse=sprintf('%04d',mouse)
    end

    fn = file_names(mouse,sess,rec)
    fileNames=[fn.rd_rec_bn '*_beh.h5'];
    filesList=dir(fullfile(fn.fold_rd_sess,fileNames));

    file=fullfile(fn.fold_rd_sess,filesList(1).name);
    %assuming only one file; read the file.
    
    dispname = strcat('reading: ',file);
    disp(dispname);
    
    
    info = h5info(file);
    n_tr = numel(info.Groups);
    
    
    table_name = info.Datasets.Name;
    table = h5read(file,['/', table_name]);
    table.odor = table.odor';

    [stream.sniff, event.breaks]  = build_continious_sniff();
    %[event.lick1,event.lick2]    = build_continious_lick();
    %[event.sniffLogic]           = sniff_logical();
    
    trPid = build_trial_structure(var, brEvents, streams);
    
    %tr = rt(tr); %runs reaction time subroutine, returns the structure.
    
%     savefilename = fullfile(fold,file);
%     savestrng = strcat('saving: ', savefilename);
%     disp(savestrng)
%     save(savefilename, 'tr');
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function [sn, breaks] = build_continious_sniff()
        % read continious sniff waveform and check the packets interuption
        % sn - continious waveform. if there is a missing packet sn is
        % performs 1st order butterworth filter
        % padded by zeros
        % breaks - [2,nb] - sequence of miisng chnnks of data onset, offset
        trEndField=var(nTrEndField).old;
        
        adc_const = 5000.0/2047; % convert factor for ADC to mV from the sniff amplifier. Only necessary for files after 20130805.
        %adc_const = 1
        
        lasTrial=find((table.(trEndField)>0),1,'last');
        preSize =table.(trEndField)(lasTrial)+1E6;
        
        disp('reading all pid packets...')
        sample_on = zeros(1,preSize);
        sn_unfiltered = nan(1,preSize);
        
        
     
        Fs = 1000;  % Sampling Frequency
     
        
        for kt = 1:n_tr
            group_name = info.Groups(kt).Name; 
            Events        = h5read(file, [group_name,'/Events']);
            sniffPackets  = h5read(file, [group_name,'/sniff']);
            
            for ie = 1:numel(Events.packet_sent_time)
                nt  = int32(numel(sniffPackets{ie}));
                ind = Events.packet_sent_time(ie) - nt-1 +(1:nt);
                if ind == 0
                    continue
                end
                sn_unfiltered(ind) = sniffPackets{ie};
                sample_on(ind) = 1;
            end
        end
        
        ind_end = ind(end)+0;
        disp (ind_end)
        disp (length(sn_unfiltered))
        sn_unfiltered = sn_unfiltered(1:ind_end);
        
        sn = zeros(1,ind_end);
        sn = sn_unfiltered*adc_const;
        
        
        breaks_on  = find(diff(sample_on)==-1);
        breaks_off = find(diff(sample_on)==1);
        
        %breaks_on = [1,breaks_on]; %starts on at first sample
        
        %checks if the first on value is greater than the first off,
        %indicating that the sample started on at 1. If the sample started off,
        %then the first on will be the first value in the array.
        if breaks_on(1) > breaks_off(1)
            breaks_on = [1,breaks_on];
        end
        
        % checks if sample is 'on' at the end, if it is, then it forces it
        % off at the last sample.
        if breaks_on(end) > breaks_off(end) 
           breaks_off = [breaks_off, ind_end];
        end
        
        nb  = numel(breaks_on);
        
%         breaks_off = breaks_off(1:nb);
        breaks = [breaks_on(:)'; breaks_off(:)'];

    end
    
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function tr = build_trial_structure(var, allEvents, streams)
        nt = numel(table.(var(1).old));
        t_pre  = 500;
        trial_record_length = 3500; %making all the trials 5 seconds in duration
        errors = [];
        trStartField=var(nTrStartField).old;   
        for kt = 1:nt,
            
            t1 = table.(trStartField)(kt) - t_pre;
            
            %checks if trialstart is a number, if not, continues.
            if ~table.(trStartField)(kt) 
                kts=num2str(kt);
                error_kt = ['error with trial ' , kts, ', skipping...'];
                disp(error_kt);
                errors=[errors,kt];
                continue
            end
                
            % record all varibles
            for kv = 1:numel(var)
                tr(kt).(var(kv).new) = table.(var(kv).old)(kt,:);
                if var(kv).relative
                    tr(kt).(var(kv).new) = tr(kt).(var(kv).new) - t1; %table.starttrial(kt);
                end
            end
            
%             tr.odor = strtrim(tr(kt).odor)
            
            % Corrections for session 1, mouse 99
            if strcmp(mouse,'ZKolfaCal') && sess==1 && tr(kt).vial==4
                tr(kt).vialConc=tr(kt).vialConc*1.;
                tr(kt).odorConc=tr(kt).odorConc*1.;
            end
            
            % Corrections for session 1, mouse 99
            if strcmp(mouse,'0099') && sess==1;
                switch tr(kt).vial
                    case 1
                        tr(kt).odor = 'Ethyl Tiglate';
                        tr(kt).vialConc=0.1;
                        tr(kt).odorConc=tr(kt).odorConc/1;
                    case 2
                        tr(kt).odor = 'Ethyl Tiglate';
                        tr(kt).vialConc=0.01;
                        tr(kt).odorConc=tr(kt).odorConc/10;
                    case 3
                        tr(kt).odor = '4-Methyl Acetophenone';
                    case 4
                        tr(kt).odor = 'Pinene';
                    case 5
                        tr(kt).odor = 'Acetophenone';
                    case 6
                        tr(kt).odor = 'Methyl Salycilate';
                        tr(kt).vialConc=0.1;
                    case 7
                        tr(kt).odor = 'Methyl Salycilate';
                        tr(kt).vialConc=0.01;
                        tr(kt).odorConc=tr(kt).odorConc/10;
                    case 8
                        tr(kt).odor = 'Eugenol';
                end
            end
            
            
            %t1 = tr(kt).starttrial - t_pre;

            t2 = tr(kt).trialstart + trial_record_length;
            
            %checks if the trial should end after the streams
            if t2>numel(stream.(streams{1}))
                kts=num2str(kt);
                error_kt = ['streams ended before trial ' , kts, ', ending trial structure building...'];
                disp(error_kt);
                errors=[errors,kt];
                tr(kt)=[];
                break
            end
            %t = t2- t1 + 1;

            
            % record all events
            for ke = 1:numel(allEvents)
                ind = (event.(allEvents{ke})(1,:)>=t1)&(event.(allEvents{ke})(2,:)<=t2);
                tr(kt).(allEvents{ke}) = event.(allEvents{ke})(:,ind) - double(t1);  %double(tr(kt).starttrial);
            end
            
            % all streams
            for ks = 1:numel(streams)
                tr(kt).(streams{ks}) = stream.(streams{ks})(t1:t2); 
            end
            
            
        end
        tr(1) = [];
    end

end

function [ avg_pid_ts,mean_timeaveraged_amplitude, mean_firstSniff_amplitude ] = pid_graph( pidStruct, odor, odorConc, vialConc, dillution, N2_pidStruct, suppress_graphs)
%PID_GRAPH analyses a single odor,concentration,vialconcentration pair.
%   [average_pid_trace,mean_amplitude] = PID_graph (input_structure, ... 
%   odor, odorConc, vialConc, suppress_graphs).
%   
%   - Input structure is derived from the 'PID_file_read.m' function.
%   - Odor is a string representing the odor for lookup
%   - odorConc is the theoretical concentration (vial concentration * air
%   dilution
%   - vialConc is the vial concentration (this is calculated based on the
%   theoretical concentration provided by the H5 record file.
%   - suppress graphs allows suppression of graph output when calling this
%   function externally.
%
%   This finds all of the trials in the input structure matching the input
%   criteria, and graphs them for quality control. It preforms DC offset correction, 
%   finds amplitude of the PID signal, discards outlier, trials 
%   (outside 2 * SD of time averaged PID amplitude), and returns an average
%   PID signal trace and an average amplitude value for the trial set.
%
%   This relies on the timing set by pid_test_protocol_0001.py. Changes to
%   these parameters requires modification to the averaging windows.
%
%   Written by CW, 29 August 2013; modified by zk Dec 03, 2013.

if nargin<6
    suppress_graphs=0;
end
len = numel(pidStruct);

% calculate the vial concentration if it doesn't exist within the structure.
% This assumes a total flow of 1000 ml/min!!!!!!

if ~isfield(pidStruct,'vialConc')
    for i = 1:len
        pidStruct(i).vialConc = pidStruct(i).liqDilution ;
    end
end

% find all of the trials meeting the input criteria.
tr_selected = (abs([pidStruct.odorConc]-odorConc)<1E-5) & (abs([pidStruct.vialConc]-vialConc)<1E-5) & strncmp({pidStruct.odor},odor,length(odor)) & [pidStruct.dillution]==dillution ;
tr_index = find(tr_selected);

%remove the trials with bad streams
badTrials = arrayfun(@(x) sum([x.sniff]==0)>3,pidStruct(tr_index));
tr_index(badTrials)=[];


if sum(tr_index)<1
    avg_pid_ts=nan;
    mean_timeaveraged_amplitude=nan;
    mean_firstSniff_amplitude=nan;
    warning('No good trials found to satisfy conditions: odor %s, [theo] = %0.1d, from [liquid] = %0.1d, extra dil 1:%d', ...
        strtrim(odor{1}),odorConc,vialConc,dillution)
    return
end

%if no trial fulfills the condition, return NaN empty avg_pid_ts,
%mean_timeaveraged_amplitued


%%%%%%%%%
% Compile all of the odor trials' sniff recordings into a matrix.
numtrials = numel(tr_index);
pid_traces = zeros(length(pidStruct(1).sniff),numtrials);

for kt = 1:numtrials
    index = tr_index(kt);
    pid_traces(:,kt) = double(pidStruct(index).sniff);
end

% DC offset correction by subtracting average of pre and post odor, plot traces
if (odorConc-0.0051<1e-5)
    fprintf('its this')
end

offset = mean(pid_traces([1:pidStruct(tr_index(1)).fvOnTime],:));
%subtract through all the presentation a straight line with the slope of the
%baseline prior to onset of fv.

baseLineTimes = [1:pidStruct(tr_index(1)).fvOnTime]';
for kt = 1:numtrials
    baseLineFit = polyfit(double(baseLineTimes),pid_traces(baseLineTimes,kt),1);
    baseLineExtended = polyval(baseLineFit,1:length(pid_traces(:,kt)));
    pid_traces(:,kt) = pid_traces(:,kt) - baseLineExtended';
end


%%%%%%%%%
% Repeat procedure for N2 signal recording
if ~isempty(N2_pidStruct)

    N2_tr_index = 1:numel(N2_pidStruct);
    %remove the trials with bad streams
    badTrials = arrayfun(@(x) sum([x.sniff]==0)>3,N2_pidStruct);
    N2_tr_index(badTrials)=[];
    
    numtrials = numel(N2_tr_index);
    N2_pid_traces = zeros(length(N2_pidStruct(1).sniff),numtrials);
    
    for kt = 1:numtrials
        index = N2_tr_index(kt);
        N2_pid_traces(:,kt) = double(N2_pidStruct(index).sniff);
    end
    
    offset = mean(N2_pid_traces([1:N2_pidStruct(N2_tr_index(1)).fvOnTime],:));
    %subtract through all the presentation a straight line with the slope of the
    %baseline prior to onset of fv.
    
    baseLineTimes = [1:N2_pidStruct(N2_tr_index(1)).fvOnTime]';
    for kt = 1:numtrials
        baseLineFit = polyfit(double(baseLineTimes),N2_pid_traces(baseLineTimes,kt),1);
        baseLineExtended = polyval(baseLineFit,1:length(N2_pid_traces(:,kt)));
        N2_pid_traces(:,kt) = N2_pid_traces(:,kt) - baseLineExtended';
    end
    average_N2_PIDtrace = mean(N2_pid_traces,2);
end


if ~suppress_graphs
    close all
    figure
    plot(pid_traces)
    
    titlestr = sprintf('PID for %s, [theo] = %0.1d, from [liquid] = %0.1d, extra dil 1:%d', ...
        odor{1},odorConc,vialConc,dillution);
    
    title(titlestr);
    ylabel('mV')
    xlabel('ms');
    
    figure
    
    surf(pid_traces,'LineStyle','none');
    title(titlestr);
    ylabel('mV')
    xlabel('ms');
    zlabel('Trial number');
    % x_on_line = [2500,4000];
    % y_on_line = [0,0];
    % annotation('line',x_on_line,y_on_line)
    
end

% find the amplitude of the PID signal, find outliers and remove them.


%the segment for average is from 500 after onset of fv to 1000 before
%closing OR last saved value
endSegment = pidStruct(tr_index(1)).fvOnTime + pidStruct(tr_index(1)).fvdur - 1000;
measureSegment=pidStruct(tr_index(1)).fvOnTime + (500 : min(endSegment,size(pid_traces,1)));

if ~isempty(N2_pidStruct)
    pid_traces = pid_traces - average_N2_PIDtrace*ones(1,size(pid_traces,2));
end

amplitude = mean(pid_traces(measureSegment,:));  % assumes that the baseline is 0
deviation = std(amplitude);
pid_traces = pid_traces(:,(abs(amplitude- mean(amplitude))<deviation*2));  % remove outlier trials
amplitude = amplitude(abs(amplitude- mean(amplitude))<deviation*2);
deviation = std(amplitude);

mean_timeaveraged_amplitude = mean(amplitude);

mean_firstSniff_amplitude = mean(mean(pid_traces(pidStruct(tr_index(1)).fvOnTime+150:pidStruct(tr_index(1)).fvOnTime+550,:)));

pid_ts = pid_traces' ;

avg_pid_ts= mean(pid_ts,1);
std_ts  = std(pid_ts);
high_pid_ts = avg_pid_ts + abs(std_ts);
low_pid_ts = avg_pid_ts - abs(std_ts);


if ~suppress_graphs
    figure;
    hold on;
    plot(high_pid_ts,'Color',[.85,0.85,0.85]);
    plot(low_pid_ts, 'Color',[.85,0.85,0.85]);
    plot (avg_pid_ts,'Color','k');
    title(titlestr);
    ylabel('mV');
    xlabel('ms');
end


% find the time-averaged deviation of each trial from the overall average.
% Make a scatter plot of these deviations to determine if there is a
% long-term shift

if ~suppress_graphs
    figure
    
    plot (amplitude)
    
    ylabel('mV')
    xlabel('trial')
    title(titlestr);
    
    top = max(amplitude)*1.1;
    last = length(amplitude);
    
    axis([1,last,0,top]);
    
    % what I want from the end is a matrix with all the samples and all the
    % trials meeting these criteria. These matrices can be graphed as line
    % plots
    
end
end

function [concPlot]=pid_multi_graph (pidStruct,odor,odorConc,vialConc,dillution, trPidN2 )
%PID_multi_graph allows comparison of multiple odor concentrations
%   
%   PID_multi_graph ( inputstruct,odor,odorConc,vialConc )
%
%   - Input structure is derived from the pid_get function.
%   - Odor is a string representing the odor for lookup
%   - odorConc is a matrix of the theoretical concentrations to be graphed
%   (read from the structure, which is read from the H5 file).
%   - vialConc is a matrix containing the vial concentrations to be
%   graphed. This is calculated from the reported concentration from the H5
%   file based on the (theoretical concentration) / (air dilution).
%
%   This program retrives sets of trials from the input structure with the
%   same input parameters. It collects the averages of these sets and
%   graphs them, allowing for comparison between vials and for assessment
%   of odor set concentration integrity. Sets with the same odorConc but
%   different vialConc are handled separately and graphed separately. If an
%   odorConc is not represented by an input vialConc, the combination is
%   ignored.
%
%   For example, if you want to graph an odorConc that is overlapping by 2
%   vials, you should input a scalar odorConc and a matrix of both
%   vialConc. If you want to graph several odorConc, some of which are
%   overlapping between vials, input the desired odorConc as a matrix with
%   each value represented once, and input all of the vialConc. The
%   overlapping values will be graphed twice, once for each vial, and the
%   non-overlapping values will be plotted only for the representitive
%   vials.
%
%   This program calls the PID_graph.m function. PID_graph.m is expecting
%   parameters defined by pid_test_protocol_0001.py. Changes to these
%   recording parameters requires changing the averaging windows defined in
%   PID_graph.m.
%
%   Written by CW, 29 August 2013.

PID_timeseries = [];
amplitude = NaN(length(vialConc)*length(dillution),length(odorConc));
amps = [];
vialConcVec=[];
dilVec=[];
vialVec=[];
leg = {};
leg2 = {};

for i = 1:numel(dillution)
    for j=1:length(vialConc)
        leg2{end+1} = sprintf('[vial] = %0.3d, ext. dil. 1:%d',vialConc(j),dillution(i));
        vialConcVec = [vialConcVec vialConc(j)];
        dilVec      = [dilVec dillution(i)];
        vialVec     = [vialVec vial_lookup(pidStruct,odor,vialConc(j))]
    end
end



%trials that have this odor
% vials   = unique([inputstruct(strncmp({inputstruct.odor},odor,length(odor))).vial]);
% trVials = inputstruct([inputstruct.vials]==vials);

for di=1:numel(dillution)
    for oC = 1:length(odorConc)
        for vC = 1:length(vialConc)
            PID_ts_i = nan;
            amplitude_i= nan;
            %if that concentratio exists for that vC; gets the plot; otherwise
            troCvC=(abs([pidStruct.odorConc]-odorConc(oC))<1E-5) & (abs([pidStruct.vialConc]-vialConc(vC))<1E-5) & strncmp({pidStruct.odor},odor,length(odor)) & [pidStruct.dillution]==dillution(di) ;
            sum(troCvC)
            if sum(troCvC)
                
                % Get signals from N2 traces to send with call to pid_graph function
                if ~isempty(trPidN2)
                    odorN2flow = pidStruct(find(troCvC,1)).nitrogen;
                    N2_pidStruct = trPidN2([trPidN2.nitrogen]==odorN2flow);
                else 
                    N2_pidStruct = [];
                end
                try
                    [PID_ts_i,amplitude_i] = pid_graph(pidStruct,odor,odorConc(oC),vialConc(vC),dillution(di),N2_pidStruct,1);
                catch
                    aaa=11;
                end
                
                if ~isnan(PID_ts_i(1))
                    PID_timeseries(:,end+1) = PID_ts_i;
                    amplitude((di-1)*(length(vialConc))+vC,oC) = amplitude_i;
                    if iscell(odor)
                        odorName=odor{1};
                    else
                        odorName=odor;
                    end
                    leg{end+1} = sprintf('[c] = %0.5d,[vial] = %0.3d, ext. dil 1:%d',odorConc(oC),vialConc(vC),dillution(di));
                end
            end
        end
        
    end
end

%fit straight lines to the concentration series of every vial
%(only for dillution=1)

amplitudeFit = nan(size(amplitude));
for vC = 1:length(vialConc)
    goodAmps =  find(~isnan(amplitude(vC,:)));
    cp = polyfit(odorConc(goodAmps),(amplitude(vC,goodAmps)),1);
    amplitudeFit(vC,:) = polyval(cp,(odorConc));
end

figure
plot(PID_timeseries)
legend(leg,'Location','best')
xlabel ('mV')
ylabel ('t(msec)')
title(odorName)

figure
plot((odorConc),amplitude,'o-')
legend(leg2,'Location','best')
ylabel('log mv')
xlabel('log concentration')
title(odorName)
hold on
plot((odorConc),amplitudeFit,'.-')

concPlot.leg   = leg2;
concPlot.vialConc = vialConcVec;
concPlot.dil   = dilVec;
concPlot.vial  = vialVec;
concPlot.amp   = amplitude;
concPlot.conc  = odorConc;

end

function [concPlot]=pid_max_graph (pidStruct,odorInfo,trPidN2)
%PID_max_graph plots the higher concentration of each vial (no dillution)
%   
%   PID_multi_graph ( inputstruct,odor,odorsInfo)
%
%   - Input structure is derived from the pid_get function.
%   - Odor is a string representing the odor for lookup
%   - odorConc is a matrix of the theoretical concentrations to be graphed
%   (read from the structure, which is read from the H5 file).
%   - vialConc is a matrix containing the vial concentrations to be
%   graphed. This is calculated from the reported concentration from the H5
%   file based on the (theoretical concentration) / (air dilution).
%
%   This program retrives sets of trials from the input structure with the
%   same input parameters. It collects the averages of these sets and
%   graphs them, allowing for comparison between vials and for assessment
%   of odor set concentration integrity. Sets with the same odorConc but
%   different vialConc are handled separately and graphed separately. If an
%   odorConc is not represented by an input vialConc, the combination is
%   ignored.
%
%   For example, if you want to graph an odorConc that is overlapping by 2
%   vials, you should input a scalar odorConc and a matrix of both
%   vialConc. If you want to graph several odorConc, some of which are
%   overlapping between vials, input the desired odorConc as a matrix with
%   each value represented once, and input all of the vialConc. The
%   overlapping values will be graphed twice, once for each vial, and the
%   non-overlapping values will be plotted only for the representitive
%   vials.
%
%   This program calls the PID_graph.m function. PID_graph.m is expecting
%   parameters defined by pid_test_protocol_0001.py. Changes to these
%   recording parameters requires changing the averaging windows defined in
%   PID_graph.m.
%


%select the max concentration for each vial
dillution = 1;
vialConc = odorInfo.vialConc;
%get the max concentration for each vial
vials = unique([odorInfo.trPid.vial]);
odorConc = arrayfun(@(x) max([odorInfo.trPid([odorInfo.trPid.vial]==x).odorConc]),vials);
vialConc = arrayfun(@(x) unique([odorInfo.trPid([odorInfo.trPid.vial]==x).vialConc]),vials);
maxConc = max(odorConc);
odor = odorInfo.odor;

PID_timeseries = [];

amplitude = NaN(length(vialConc)*length(dillution),length(odorConc));
amp_snf1 =  NaN(length(vialConc)*length(dillution),length(odorConc));
vialConcVec=[];
dilVec=[];
vialVec=[];
ratioVec = [];

leg = {};
leg2 = {};

for i = 1:numel(dillution)
    for j=1:length(vialConc)
        leg2{end+1} = sprintf('[vial] = %0.3d, ext. dil. 1:%d',vialConc(j),dillution(i));
        vialConcVec = [vialConcVec vialConc(j)];
        dilVec      = [dilVec dillution(i)];
        vialVec     = [vialVec vial_lookup(pidStruct,odor,vialConc(j))]
    end
end



%trials that have this odor
% vials   = unique([inputstruct(strncmp({inputstruct.odor},odor,length(odor))).vial]);
% trVials = inputstruct([inputstruct.vials]==vials);

for di=1:numel(dillution)
    for oC = 1:length(odorConc)
        for vC = 1:length(vialConc)
            PID_ts_i = nan;
            amplitude_i=nan;
            amp_snf1_i=nan;
            %if that concentratio exists for that vC; gets the plot; otherwise
            troCvC=(abs([pidStruct.odorConc]-odorConc(oC))<1E-5) & (abs([pidStruct.vialConc]-vialConc(vC))<1E-5) & strncmp({pidStruct.odor},odor,length(odor)) & [pidStruct.dillution]==dillution(di) ;
            sum(troCvC)
            
            
            if sum(troCvC)
                ratio = maxConc/odorConc(oC);
                ratioVec = [ratioVec ratio];
                
                % Get signals from N2 traces to send with call to pid_graph function
                if ~isempty(trPidN2)
                    odorN2flow = pidStruct(find(troCvC,1)).nitrogen;
                    N2_pidStruct = trPidN2([trPidN2.nitrogen]==odorN2flow);
                else
                    N2_pidStruct = [];
                end
                
                [PID_ts_i,amplitude_i,amp_snf1_i] = pid_graph(pidStruct,odor,odorConc(oC),vialConc(vC),dillution(di),N2_pidStruct,1);
            end
            
            
            if ~isnan(PID_ts_i)
                PID_timeseries(:,end+1) = PID_ts_i*ratio;
                amplitude((di-1)*(length(vialConc))+vC,oC) = amplitude_i;
                amp_snf1((di-1)*(length(vialConc))+vC,oC) = amp_snf1_i;
                if iscell(odor)
                    odorName=odor{1};
                else
                    odorName=odor;
                end
                leg{end+1} = sprintf('[c] = %0.5d,[vial] = %0.3d, ext. dil 1:%d',odorConc(oC),vialConc(vC),dillution(di));
            end
            
        end
    end
end

figure
plot(PID_timeseries,'.')
legend(leg,'Location','best')
xlabel ('mV')
ylabel ('t(msec)')
title(odorName)

concPlot.leg   = leg2;
concPlot.vialConc = vialConcVec;
concPlot.dil   = dilVec;
concPlot.vial  = vialVec;
concPlot.amp   = amplitude;
concPlot.amp_snf1 = amp_snf1;
concPlot.conc  = odorConc;

end

function [odors]=pid_all_series(trPid)
%make a structure with all the sets of stimuli,
%and plot and save the matrix of amp of response vs. concentration
%odor
 % all the concentrations, all the liquid concentrations, all the
 % dillutions
 
 %list of odors
 odorList=unique({trPid.odor})
 
 %all the trials for an odor
 for io=1:numel(odorList)
     if strncmp('one',odorList(io),4)
         continue
     end     
     if any(strncmpi('none',odorList,4))
         trialsN2 = find(strcmp(odorList(strncmpi('none',odorList,4)),{trPid.odor}));
         trialsN2 = trialsN2;
         trPidN2  = trPid(trialsN2);
     else
         trPidN2 = [];
     end
     
     trials=find(strcmp(odorList(io),{trPid.odor}));
     
     odors(io).odor     = odorList(io)
     odors(io).trials   = trials;
     odors(io).trPid    = trPid(trials);
     odors(io).odorConc = unique([trPid(trials).odorConc]);
     odors(io).vialConc = unique([trPid(trials).vialConc]);
     odors(io).dillution= unique([trPid(trials).dillution]);
     odors(io).concPlot = pid_multi_graph( trPid, odors(io).odor, odors(io).odorConc, odors(io).vialConc, odors(io).dillution, trPidN2);
     odors(io).maxPlot  = pid_max_graph( trPid, odors(io), trPidN2);
     %if there are more than 1 viales with the same liquid dillution,
     %calculate their effective dillution (relative to hightest conc).
     %WORKS WITH ONLY TWO MATCHED VIALS!!!
     odors(io).equalConc= find(diff([odors(io).concPlot.amp(1:numel(odors(io).concPlot.vialConc),:)])<1E-6);
     concRate=ones(1,numel(odors(io).vialConc));
     if ~isempty(odors(io).equalConc)
         [sortedC,indC]=sort(odors(io).concPlot.vialConc);
         concRate(1:end-1)=arrayfun(@(x) odors(io).concPlot.amp(x,odors(io).equalConc)/odors(io).concPlot.amp(indC(end),odors(io).equalConc),indC(1:end-1))
     end
     odors(io).concRate=concRate;
 end
 
end

function [ outstruct ] = pid_graph_all( inputstruct )
% Groups trials by odor & concentration and plots.
%   This allows a quick check for drift and 

outstruct = struct;

%remove all the traces that have breaks
noBreaks=arrayfun(@(x) isempty(x.breaks),inputstruct);
inputstruct(~noBreaks)=[];
outstruct.Trials = inputstruct;

len = length(inputstruct);
outstruct.conc = [];
outstruct.od = {};



%%make a vector for the odors and concentrations to use for masking trials
%%based on their properties.
for kt = 1:len
    
    outstruct.od(end+1)   = cellstr(inputstruct(kt).odor);
    outstruct.conc(end+1)  = inputstruct(kt).odorConc;

end

outstruct.concentrationlist = unique(outstruct.conc);
outstruct.odorlist = unique(outstruct.od);


% make indeces of each odor:concentration pair. Make a structure element f

iO = [];
iC = [];



for ko = 1:length(outstruct.odorlist)
    
    odorvar = genvarname(outstruct.odorlist{ko});
    outstruct.(odorvar) = {};
    
    
    for kt = 1:len
        if isequal(outstruct.odorlist(ko),outstruct.od(kt))
            iO(kt) = 1;
        else
            iO(kt) = 0; 
        end 
    end
    
    iO = logical(iO);
    
    for kc = 1:length(outstruct.concentrationlist)
        
        iC = (outstruct.conc == outstruct.concentrationlist(kc));
        
        intersection = iO .* iC;
        outstruct.(odorvar){kc} = find(intersection);
   
    end
end


% great. now you've done that easy stuff. Let's plot each individual trace, making a figure for each odor:concentration pair. This will help us determine if the odor signal was consient through the session.

for ko = 1:length(outstruct.odorlist)
    
    odorvar = genvarname(outstruct.odorlist{ko});
    
    for kc = 1:length(outstruct.(odorvar))
        
        ind = outstruct.(odorvar){kc};
        if isempty(outstruct.(odorvar){kc})
            continue
        else
            sniffmat = [];
            
            
            for kts = 1:length(ind)
                %remove the traces that have lost packets
                rawTrace = outstruct.Trials(ind(kts)).sniff;
                if any(isnan(rawTrace))
                    rawTrace = nan(size(rawTrace));
                end
                sniffmat(:,kts) = rawTrace;
                
            end
            
            
            
            figure
            lentrial = length(sniffmat);
            x = linspace(1,lentrial,lentrial);
            offset = mean(sniffmat([1:500],:));
            
            for offf = 1:length(offset)
                sniffmat(:,offf) = sniffmat(:,offf) - offset(offf);
            end
            
            plot (x,sniffmat)
            titlestr1 = sprintf(' %d concentration ',outstruct.concentrationlist(kc));
            titlestr2 = [odorvar, titlestr1];
            title(titlestr2)
            
            figure
            surf(sniffmat,'LineStyle','none');
            ylabel('mV')
            xlabel('ms');
            zlabel('Trial number');
            title(titlestr2)
        end
    end
    
end
end

function [pidAn]=pid_analysis(pid,sess,stat)
%pulls pid data from a station
%gets the pid (pid_get)
%runs pid_all_series to view the signals

%pid files end with _pid.h5: make sure voyeur records like that

%pull the data if stat was entered
if nargin>2 && ~isempty(stat) && ~strcmp(stat,'local')
    pull_data(pid,sess,stat);
end

%get the pid. make it always rec_a, don't complicate it
%(otherwies get filenames, dir, count the recs, and go rec by rec.
recList=['a'];
for ir=1:numel(recList)
    trPid=pid_get(pid,sess,recList(ir));
    pidAn(ir).rec=recList(ir);
    pidAn(ir).sess=sess;
    pidAn(ir).odors=pid_all_series(trPid);
end


% This will do the little calculation of concentration corrections for you
%   Note: assumes that i) only 2 vials compared, and ii) original flow rates were 100:900
% for io = 1:numel(pidAn.odors)
%     concRatio = pidAn.odors(io).maxPlot.conc(2)/pidAn.odors(io).maxPlot.conc(1);
%     amp = pidAn.odors(io).maxPlot.amp_snf1;
%     if concRatio<1
%         amp = [amp(:,2) amp(:,1)];
%         amp = [amp(2,:); amp(1,:)];
%         concRatio = pidAn.odors(io).maxPlot.conc(1)/pidAn.odors(io).maxPlot.conc(2);
%     end
%     ampRatio = (amp(1,1)*concRatio)/amp(2,2);
%     newAirFlow = round(900*ampRatio);
%     newN2Flow  = 1000-newAirFlow;
%     if abs(ampRatio-1)<0.2 && newN2Flow<100 && newN2Flow>10
%         fprintf('\n%s\n  (10*vial_low)/vial_high ratio around time of first sniff is %4.4f \n  Correct lower conc vial flows to %i:%i\n',strtrim(pidAn.odors(io).odor{1}),ampRatio,newN2Flow,newAirFlow)
%     else
%         newVialConc = pidAn.odors(io).maxPlot.vialConc(1)/ampRatio;
%         fprintf('\n%s\n  (10*vial_low)/vial_high ratio around time of first sniff is %4.4f \n  CHANGE LIQUID DILUTIONS!\n  Try new lower vial conc %4.4f (current is %4.4f)\n',strtrim(pidAn.odors(io).odor{1}),ampRatio,newVialConc,pidAn.odors(io).maxPlot.vialConc(1))
%     end
% end
end





function vial = vial_lookup(trPid,odor,vialConc)
%lookup the vial that has that odor at that concentration, to 6 digits
trialsThisOC=find(round(100000*[trPid.vialConc])==round(100000*vialConc) & strncmp({trPid.odor},odor,length(odor)));
vial=unique([trPid(trialsThisOC).vial]);

if numel(vial)>1
    warning('There were %d vials of %s with concentration %f \n', numel(vial), odor, vialConc);
    vial=vial(1);
end

if isempty(vial)
        error('There were no vials of %s with concentration %f \n', odor, vialConc);
end

end

function [gf, gs] = subplot_fig(ny, nx, fig)
    if exist('fig', 'var')
        gf = figure(fig); clf
    else
        gf = figure;
    end
    
    dx = 0.9/nx;    dy = 0.9/ny;
    wx = dx*0.95;   wy = dy*0.95;
    for ix = 1:nx
        for iy = 1:ny
            gs(iy,ix) = subplot('position', [0.08+(ix-1)*dx, 0.95-iy*dy, wx, wy]);
        end
    end
end
