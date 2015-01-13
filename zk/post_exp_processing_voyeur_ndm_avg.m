% Data preparation for spike sorting
% Dima Rinberg, August 2013
%
% ZK - Jan 2014:
% DO NOT LEAVE THIS SCRIPT BROKEN.
% IT IS CALLED BY THE PREPROCESSING SCRIPTS IN THE COMPUTATIONAL SERVER
% IF YOU ARE GOING TO PLAY WITH IT, MAKE ANOTHER BRANCH.
% SEE DETAILED DESCRIPTION IN EVERNOTE.
%
%2013/09/12
%This is the latest debugged, working version that does not push data to
%server.
% Version 2013/09/09, modifications to v04
% function (and function fn) now takes 4 parameters:
% stat is for the station where the data is stored.
% if stat is not input, local source is assumed and no data is pulled.
% function fn will not generate fnames for rec if input in rec is empty
% ('').
% added function pull_data(mouse, sess, stat) to pull a session from the location
% where a station. is mounted.
%in line 171, 
% do not inherit run structure fields from prev. fields (in function
% function [info, rec] = read_log_file())
% 2013/12/09
% This is the program that was in
% ../preprocessingscripts/post_exp_preprocessing_voyeur_working
% Added field group to chan structures, which is read from the new
% channe_config_file.txt with 4 columns.
% Modified (corrected) read_chan_config so that numbers and names of
% channels are not matched by index in the file but by number declared in
% the config_file.

function post_exp_processing_voyeur_ndm_avg(mouse, sess,stat)
global pp
%add the package of common files and folders to the path.
%here are the functions for file_naming, for instance, and any other
%functions that several script bundles will use.
%that is /baseFolder/ephysDataManagement/current/include
%if it's called from the file server, the MATLABPATH already includes them,
%so no include is necessary.
[~,computerName]=system('hostname');
if ~strcmp(strtrim(computerName),'flipper')
    includePath=fullfile(fileparts(pwd),'current','include');
    addpath(includePath);
end

% Functions of include that it uses:
%   - file_names(mouse,sess,rec,stat)

pp.read_log_file  = @read_log_file;
pp.read_data_dir  = @read_data_dir;
pp.ss_prep        = @ss_prep;
pp.xml_make       = @xml_make;
pp.resampling     = @resampling;
pp.trial_prep     = @trial_prep;
pp.basicVoyerRead = @basicVoyerRead;

if ~exist('mouse', 'var')    
    mouse = 'SDbehavingM72';
    sess  = '001';
    %stat = '01';
end

if nargin>2 && ~isempty(stat) && ~strcmp(stat,'local')
    pull_data(mouse,sess,stat);
end

%read_info(mouse, sess);
%ss_prep(mouse, sess)
xml_make(mouse,sess);
%resampling(mouse, sess)
%trial_prep(mouse, sess)
end

function read_info(mouse, sess)
% the program read the following inofmration:
% 1 - log file
% 2 - raw data file directory
% 3 - electrode channel config
% 4 - meta data for individual recordings
% It prints the structure rec for debugging.
% It saves the data into ss_fold, sess_info.mat file
global info;

    fn = file_names(mouse, sess)

    [info, rec] = read_log_file();
 
    rec   = read_data_dir(rec);
    rec   = read_meta_data(rec);

    ff = fieldnames(rec);
    id = find(strcmp(ff, 'chan'), 1);
    ff = ff([1:id-1,id+1:end, id]);
    id = find(strcmp(ff, 'run'), 1);
    ff = ff([1:id-1,id+1:end, id]);
    rec = orderfields(rec, ff);

    if ~exist(fn.fold_ss_sess, 'dir')
        mkdir(fn.fold_ss_sess)
    end
    info.rec = rec;
    save(fn.ss_sess_info, 'info')

    for irec=1:numel(rec)
        disp(['rec' num2str(irec)]);
        disp('-----');
        disp(rec(irec));
        for irun=1:numel(rec(irec).run)
            disp(['run' num2str(irun)]);
            disp('------');
            disp(rec(irec).run(irun));
        end
        disp('*--------------------------------------*');
    end
    
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function [info, rec] = read_log_file()
        info   = struct();
        rec    = struct();
        irec   = 0;
        if ~exist(fn.log, 'file')
            error('no log file')
        end
        fid = fopen(fn.log, 'r');

        neof   = true;
        plevel = 1;
        while neof
            % read a line from log file
            tl = fgetl(fid);
            % if end of the file, finish the prorgam, 
            % if the line starts with ':' - read the parameter
            if tl== -1
                neof =0;
            elseif ~isempty(tl)&&(tl(1)==':')
                % extract parameter name 'var' and its value 'value'
                [var, value] = strtok(tl(2:end), ':');
                var   = strtrim(var);
                value = strtrim(strtok(value(2:end), '%'));
                nvalue = str2double(value);

                
                switch var
                    case 'rec'
                        plevel = 2;
                        irec   = irec + 1;
                        if irec > 1
                            rec(irec) = rec(irec-1);
                        end
                        rec(irec).name = value
                        % do not inherit run structure fields from prev. fields
                        if isfield(rec(irec),'run') && ~isempty(rec(irec).run)
                            rec(irec).run=[];
                        end
                            
                    case 'run'
                        plevel = 3;
                        irun = nvalue;
                        if isnan(irun)
                            error('wrong run number')
                        end
                        rec(irec).run(irun).num = irun;
                     otherwise
                        % record parameter value, if it is a numerical convert to numerical
                        if ~isnan(nvalue)
                            value = nvalue;
                        end
                        
                        switch plevel
                            case 1
                                info.(var) = value;
                            case 2
                                rec(irec).(var) = value;
                            case 3
                                rec(irec).run(irun).(var) = value;
                        end
                end
            end
        end
        fclose(fid);
%         save('recdbg.mat','rec')
    end
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rec = read_data_dir(rec)
        % the program reads raw data sess folder, sort files by records and runs,
        % and reads meta information for records and runs
        % the program uses the 'rec' structure saved in 'ss_data\se_xxx_xx' folder
        % it uses the paprmaters from these structure
        % if it finds a record with a new name it uses the paarmeters of the last
        % record.
        % the output is saved back to the same file, same structure

        ftype  = {'ephys_data', 'ephys_meta', 'behav_data'};

        % total number of runs for each rec
        nr = zeros(1,numel(rec));

        for kf = 1:numel(ftype)
            % read the list of file of the same type 'ftype' from raw data session
            % folder
            fnam = fullfile(fn.fold_rd_sess, ['*',info.(ftype{kf})]);
            ll = dir(fnam);

            for il = 1:numel(ll)
                % defining the record
                [rname, rem] = strtok(ll(il).name, '_');
                % finding if the name of the rec exists, or create a new
                % rec if not
                irec = find(strcmp(rname, {rec.name}), 1);
                %disp(['Rec ' num2str(irec)])
                if isempty(irec)
                    % record does not exist
                    irec           = numel(rec)+1;
                    rec(irec)      = rec(end);
                    rec(irec).name = rname;
                    nr(irec)       = 0;
                    warnglg('new record name')
                end

                % defining the run number
                run = str2double(strtok(rem(2:end),'_'));
               % disp(['nr(' num2str(irec) ')= ' num2str(nr(irec))]);
                if ~isnan(run)
                    if nr(irec) == 0
                        irun     = 1;
                    else
                        lastIfound=irun;
                        irun = find([rec(irec).run.num]==run,1);
                        if isempty(irun)
                            irun = lastIfound+1;
                        end
                    end
                    rec(irec).run(irun).num = run;
                    rec(irec).run(irun).(ftype{kf}) = ll(il).name;
                    nr(irec) = length([rec(irec).run.num]);
                end
            end
        end
    end

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rec = read_meta_data(rec)
        for kr = 1:numel(rec)

            fn_meta = fullfile(fn.fold_rd_sess, rec(kr).run(1).ephys_meta);
            meta    = read_meta_data(fn_meta);

            rec(kr).nChan_rd = meta.nChans;
            rec(kr).gain     = meta.auxGain;
            rec(kr).rangeMax = meta.rangeMax;
            rec(kr).rangeMin = meta.rangeMin;

            % read channel config file
            fn_config = fullfile(fn.fold_rd_sess, rec(kr).channel_config);
            if ~exist(fn_config)
                fn_config = fullfile(fn.fold_config, rec(kr).channel_config);
            end
            chan_rd  = read_chan_config(fn_config);

            % analog channels
            in_analog = find(strcmp({chan_rd.type}, 'A'));
            nChan     = numel(in_analog);
            for kc = 1:nChan
                ic = in_analog(kc);
                chan(kc).num_rd = chan_rd(ic).num;
                chan(kc).name   = chan_rd(ic).name;
                if isfield(chan_rd,'group')
                    chan(kc).group  = chan_rd(ic).group;
                end
            end
            rec(kr).nChan = nChan;
            rec(kr).chan  = chan;

            % digital channels
            in_digital = find(strcmp({chan_rd.type}, 'D')&~strcmp({chan_rd.name}, 'PLtrig'));
            nDigChan   = numel(in_digital);
            if nDigChan > 0
                for kc = 1:nDigChan
                    ic = in_digital(kc);
                    digChan(kc).num_rd = chan_rd(ic).num;
                    digChan(kc).name   = chan_rd(ic).name;
                end
                rec(kr).digChan  = digChan;
            end
            rec(kr).nDigChan = nDigChan;

            % PL trigger
            iPltrig=find(strcmp({chan_rd.name}, 'PLtrig'));
            rec(kr).pltrig_num_rd  = chan_rd(iPltrig).num;

            offset = 0;
            for ir = 1:numel(rec(kr).run),  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                fn_meta = fullfile(fn.fold_rd_sess, rec(kr).run(ir).ephys_meta)
                meta   = read_meta_data(fn_meta);
                nSampl = meta.fileSizeBytes/meta.nChans/2;

                rec(kr).run(ir).nSampl = nSampl;
                rec(kr).run(ir).offset = offset;
                offset = offset + nSampl;

            end,    % ir ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        end,    % kr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
            
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function q = read_meta_data(fnam)
            fid = fopen(fnam, 'r');
            if fid==-1
                error('no meta file')
            end

            neof = 1;
            while neof
                tl = fgetl(fid);
                if tl==-1
                    neof = 0;
                elseif ~isempty(strfind(tl, '='))
                    [name, rem] = strtok(tl, '=');
                    name   = strtrim(name);
                    value  = strtrim(rem(2:end));
                    nvalue = str2double(value);
                    if ~isnan(nvalue)
                        value = nvalue;
                        %fprintf('Par name %16s, value %f\n',name, nvalue);
                    end
                   q.(name) = value;
                                  
                end
            end
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function ch = read_chan_config(fnam)
            disp(['reading channel config file ' fnam '...']);
            if ~exist(fnam, 'file')
                error('no channel_config file')
            end
            fid  = fopen(fnam, 'r');
            str = textscan(fid, '%d  %s   %s   %s');
            num = double(str{1});
            type = str{2};
            name = str{3};
            numel(name);
            group= str{4};
            numel(group);
            fclose(fid);

            ch = struct();
            for kc = 1:numel(num)
                ch(kc).num  = num(kc);
                ch(kc).name = name{kc};
                ch(kc).type = type{kc};
                %backwards compatibility: if there is no group column, the
                %first element of group is empty. then group field should
                %not exist.
                if ~isempty(group{1});
                    ch(kc).group= group{kc};
                end
            end
            disp(['done, file has description of ' num2str(numel(str{1})) ' channels.']);
        end

    end
  
end

function ss_prep(mouse, sess)
    tic
    fn = file_names(mouse, sess);
    if ~exist(fn.ss_sess_info, 'file')
        error('no sess_info file')
    end
    
    q = load(fn.ss_sess_info);
    rec = q.info.rec;

    for kr = 1:numel(rec),   % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fprintf('==================================================\n')
        fprintf('rec: %s \n', rec(kr).name)
        
        ipl      = rec(kr).pltrig_num_rd;         % PL trigegr channel number
        chan     = rec(kr).chan;                  % channel parameters
        nChan    = rec(kr).nChan;                 % number of channels
        nChan_rd = rec(kr).nChan_rd;              % number of raw data channels
        
        pl_extr = strcmp(rec(kr).pl_trigger, 'yes')&~isempty(ipl);
        
        % create a record SS folder and SS file
        fn = file_names(mouse, sess, rec(kr).name);
        mkdir(fn.fold_ss_rec)
        fid_rec = fopen(fn.ss_rec, 'w');
        
        for ir = 1:numel(rec(kr).run)
            fprintf('--------------------------------------------------\n')
            fprintf('run: %d ', rec(kr).run(ir).num)
            
            run    = rec(kr).run(ir);                  % run parameters
            
            % open raw data file
            fn_data = fullfile(fn.fold_rd_sess, run.ephys_data);
            fid_r = fopen(fn_data, 'r');
            fprintf('file: %s\n ', fn_data)
            if pl_extr
                fprintf('pl_trigger \n')
                pl   = read_pl_trigger(fid_r, ipl, nChan_rd);
                if isempty(pl)
                    pl_extr = 0;
                else
                    pl_freq = rec(kr).sampling_freq/mean(diff(pl));
                    if (pl_freq < 50)||(pl_freq > 70)
                        pl_extr = 0;
                    end
                end
            end
            
            for kc = 1:nChan,
                fprintf('chan: %4d   %3d  %7s  ', kc, chan(kc).num_rd, chan(kc).name)
                X = read_analog_channel(fid_r, chan(kc).num_rd, nChan_rd);
                X = X - mean(X);
                
                if pl_extr
                    fprintf('-> pl trigger     ')
                    X = remove_pl_trigger(X, pl);
                end
                
                fprintf('-> writing data')
                write_analog_signal(X, fid_rec, run.offset, kc, nChan)
                fprintf('   %d\n', round(toc))
            end
            fclose(fid_r);
        end
        fclose(fid_rec);
    end
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function pl = read_pl_trigger(fid, ich, nch)
        S = read_analog_channel(fid, ich, nch);
        thr = mean(S);
        pl = find(diff(S>thr)==1);
        figure
        plot(S(1:10000))
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function Y = read_analog_channel(fid, ich, nch)
        fseek(fid, 2*(ich-1), 'bof');
        Y = fread(fid, inf, 'int16', (nch-1)*2);
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function X = remove_pl_trigger(X, pl)
        if ~isempty(pl)
            %remove pl trigger
            n_templ  = ceil(mean(diff(pl))*1.2);
            n_pl     = numel(pl)-2;
            in_templ = (1:n_templ)'*ones(1,n_pl) + ones(n_templ,1)*pl(1:n_pl)'-1;
            templ    = mean(X(in_templ),2);

            for it = 1:n_pl;
                in = pl(it)+1:pl(it+1);
                nn = length(in);
                X(in) = X(in) - templ(1:nn);
            end
        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function write_analog_signal(X, fid, offset, kch, nch)
        fseek(fid, (offset*nch+kch-1)*2, 'bof');
        fwrite(fid, X(1), 'int16');
        fwrite(fid, X(2:end), 'int16', (nch-1)*2);
    end
end

function xml_make(mouse, sess,rec)
%this function reads trhough the info of the session and writes an XML file
%(for each rec) to send to ndmanager trhough the bash script.
%it reads the info file that is created in the ss_data folder, and then
%makes the xml file in the corresponding folder.

%the structure of the xml file is:
%keys either have values, or indicate sections. They also have options.
%most of the values are read from the info file; some are input here
%(as of jan 2014, hardcoded; since they don't change often).

%the names of the fields of the structure are the names of the keys in the
%file, and they follow the same tree organization.
% To indicate that a key is for a value, it has a 'val' field with the
% value.

%add the package of common files and folders to the path.
%here are the functions for file_naming, for instance, and any other
%functions that several script bundles will use.
%that is /baseFolder/ephysDataManagement/current/include
[~,computerName]=system('hostname');
if ~strcmp(strtrim(computerName),'flipper')
    includePath=fullfile(fileparts(pwd),'current','include');
    addpath(includePath);
end
% Functions of include that it uses:
%   - file_names(mouse,sess,rec,stat)
%   - number2color

%for debugging
global dbg;
%get the info structure
fprintf('\n-----------------------------------------------------------\n');
fprintf('Preparing sorting parameters for mouse %s, session %s...\n',num2str(mouse),num2str(sess));
fn = file_names(mouse, sess);

if ~ isfield(fn,'xml_template') || isempty(fn.xml_template)
    warning(sprintf('This machine is not configured for making ndmanager .xml parameter files; skipping this function (no xml template file)'));
    return
end

fprintf('\t session info file: %s\n',fn.ss_sess_info);

sInfo=load(fn.ss_sess_info);
info=sInfo.info;

if nargin<3 || isempty(rec)
    fprintf('\t no record specified, going for all of them...\n');
    recs=1:numel(info.rec);
else
    recs=strcmp({info.rec.name},rec);
end

for iRec=recs
    fprintf('\t rec %d (%s)\n',iRec,sInfo.info.rec(iRec).name);
    %info is sess_info
    rec=sInfo.info.rec(iRec);
    %if this rec doesn't have a group config; skip the xml file with a
    %warning
    if ~isfield(rec.chan,'group')
        warning('Rec %s does not seem to have sorting group information in the channel config file. I will not sort these spikes.',rec.name);
        continue
    end
    fn = file_names(mouse, sess,rec.name);
    parameters.opt=sprintf('creator="ndManager-1.2.1" version="1.0"');
    %general info
    gInfoFields = {'experimenters','description','notes'};
    generalInfo.date.val=strrep(info.date,'/','-');
    for iFi=1:numel(gInfoFields)
        if isfield(info,(gInfoFields{iFi}))
            generalInfo.(gInfoFields{iFi}).val=info.(gInfoFields{iFi});
        else
            generalInfo.(gInfoFields{iFi})=struct;
        end
    end
    if isempty(fields(generalInfo.experimenters))
        generalInfo.experimenters.val='Rinberg Lab';
    end
    parameters.generalInfo=generalInfo;
    %acquisitionsystem
    acquisitionSystem.nBits.val=16;
    acquisitionSystem.nChannels.val=rec.nChan;
    acquisitionSystem.samplingRate.val=rec.sampling_freq;
    acquisitionSystem.voltageRange.val=sprintf('%2d',rec.rangeMax-rec.rangeMin);
    acquisitionSystem.amplification.val=400;
    acquisitionSystem.offset.val=0;
    parameters.acquisitionSystem=acquisitionSystem;
    %fieldPotentials
    parameters.fieldPotentials.lfpSamplingRate.val=1250;
    %files
    parameters.files.file.samplingRate.val=rec.sampling_freq;
    parameters.files.file.extension.val='fil';
    %group parameters brought from info structure
    [chGroupInfo, chansInfo]=get_channel_groups(rec);
    %anatomicalDescription
    for ig=1:numel(chGroupInfo)
        for iCh=1:chGroupInfo(ig).nCh
            aDgroup(ig).channel(iCh).opt=sprintf('skip="0"');
            aDgroup(ig).channel(iCh).val=chGroupInfo(ig).chan(iCh);
        end
    end
    parameters.anatomicalDescription.channelGroups.group=aDgroup;
    %SpikeDetection
    props={'nSamples','peakSampleIndex','nFeatures'};
    vals=[32,16,3];
    aDS=aDgroup([chGroupInfo.sortable]==1);
    for ig=1:numel(aDS)
        sG(ig).channels=aDS(ig);
        for iChG=1:numel(sG(ig).channels.channel)
            sG(ig).channels.channel(iChG).opt=[];
        end
        for ip=1:numel(props)
            sG(ig).(props{ip}).val=vals(ip);
        end
    end
    parameters.spikeDetection.channelGroups.group=sG;
    %units
    parameters.units=[];
    %neuroscope
    neuroscope.opt=sprintf('version="1.2.5"');
    %miscellaneous
    miscellaneous.screenGain.val=0.2;
    miscellaneous.traceBackgroundImage.val=[];
    neuroscope.miscellaneous=miscellaneous;
    %video
    video.rotate.val=0;
    video.flip.val=0;
    video.videoImage.val=[];
    video.positionsBackground.val=0;
    neuroscope.video=video;
    %spikes
    neuroscope.spikes.nSamples.val=32;
    neuroscope.spikes.peakSampleIndex.val=16;
    %channels
    default_color=sprintf('#993399');
    for iCh=1:numel(chansInfo)
        chNum=chansInfo(iCh).number;
        chGrp=chansInfo(iCh).gId;
        ch.channelColors(iCh).channel.val=chNum;
        if exist('number2color')
            chColor=number2color(chansInfo(iCh).gId);
        else
            chColor=default_color;
        end
        ch.channelColors(iCh).color.val=chColor;
        ch.channelColors(iCh).anatomyColor.val=chColor;
        ch.channelColors(iCh).spikeColor.val=chColor;
        ch.channelOffset(iCh).channel.val=chNum;
        ch.channelOffset(iCh).defaultOffset.val=0;
    end
    neuroscope.channels=ch;
    
    parameters.neuroscope=neuroscope;
    
    %2014 January go upto here defining parameters.
    %all the rest just concatenates a fixed ndms script, since there are no
    %parameters for we to change; after closing the key 'neuroscope',
    
    strprint(parameters);
    %if the structure is parameter
    clear dbg;
    dbg=parameters;
    rec.chGroupInfo=chGroupInfo;
    newRec(iRec)=rec;
    
end %for iRec=recs
info.rec=newRec;
save(fn.ss_sess_info,'info');
fprintf('\tDone making ndm .xml file(s) for mouse %s, session %s!\n',num2str(mouse),num2str(sess));
fprintf('-----------------------------------------------------------\n');

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function [group, channels] =get_channel_groups(rec)
        %get all the channels: their number in the data file counting from zero, 
        %the group they belong to, their order in
        %the group, and if the group is sortable (anatomical vs spike
        %group)
        %this was read when the channel_config file was read, and is in the
        %sess_info structure in the ss_data.
        
        for iCh=1:numel(rec.chan)
            %channels(iCh).number=rec.chan(iCh).num_rd-1;
            channels(iCh).number=iCh-1; %channels as in the preprocessed file.
            %there is an old way of parsing the group information.
            if isempty(strfind(rec.chan(iCh).group,'_'))
                channels(iCh).group=rec.chan(iCh).group(1);
                channels(iCh).grPos=str2num(rec.chan(iCh).group(2:end));
            else
                [channels(iCh).group,gPos]=strtok(rec.chan(iCh).group,'_');
                channels(iCh).grPos=str2num(gPos(2:end));
            end
        end
        
        %now find the groups
        %figure out how many channels they have, and whether they are ephys
        %make a group structure compatible with the xml file with fields
        %name
        %nCh
        %chan
        %sortable
        groupList=sort(unique({channels.group}));
        for ig=1:numel(groupList)
            gr=groupList(ig);
            group(ig).name=gr;
            grChan=channels(strcmpi({channels.group},gr));
            %assign a group id in every channel
            for icg=[grChan.number]
                channels(find([channels.number]==icg)).gId=ig;
            end
            %sort by order in the group
            [grS,iGsort]=sort([grChan.grPos]);
            group(ig).chan=[grChan(iGsort).number];
            group(ig).nCh=numel(group(ig).chan);
            if grChan(1).grPos==0
                group(ig).sortable=0;
            else
                group(ig).sortable=1;
            end  
        end
    end
    %of function group=get_channel_groups(rec)
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function meta = strprint( stru,meta)
        % prints a ndm_structure in ndm xml format to a file.
        % when called from outside, there is no meta, so it initializes.
        % if the structure is parameters, it adds the ndm header <?xml
        % version='1.0'?> in level 0,
        % and pastes the standard, post-ephys parameters part of the xml file in
        % the end.
        %
        % checks how many fields the structure has.
        % for each of them:
        % if it has a val field, it prints it.
        % if none, it prints the null key
        if nargin<2
            level=0;
            struName=inputname(1);
            meta.fWid=fopen(fn.ss_xml,'w');
            fprintf(meta.fWid,'<?xml version=\''1.0\''?>\n');
        else
            level=meta.level+1;
            struName=meta.name;
        end
        
        %make the indent
        levelTab='';
        for it=1:level
            levelTab=[levelTab ' '];
        end
        
        %get the contents of the struct and print them
        keys=fields(stru);
        %if its a null value field, print it as a null value field
        if numel(keys)==0
            fprintf(meta.fWid,'%s<%s/>\n',levelTab,struName);
            
        else
            %if it is a non-null field, open the key
            %check if it has an option
            optStr='';
            if sum(strcmp('opt',keys))
                if ~isempty(stru.opt)
                    optStr=[' ' stru.opt];
                end
                keys(strcmp(keys,'opt'))=[];
            end
            
            fprintf(meta.fWid,'%s<%s%s>',levelTab,struName,optStr);
            %if it is a value, print the value and get ready to close
            if sum(strcmp('val',keys))
                fprintf(meta.fWid,'%s',num2str(stru.val));
                levelTab='';
            else
                fprintf(meta.fWid,'\n');
                keys(strcmp(keys,'val'))=[];
                %if it has other fields, print them and print carriage return
                for ik=1:numel(keys)
                    %fprintf('printkey %s<%s>\n',levelTab,keys{ik});
                    key=keys{ik};
                    subStruArray=stru.(key);
                    %             subStruArray(1)
                    for ia=1:length(subStruArray)
                        meta.name=keys{ik};
                        meta.level=level;
                        subStru=subStruArray(ia);
                        meta=strprint(subStru,meta);
                    end
                    
                end
            end
            if ~strcmp(struName,'parameters')
                fprintf(meta.fWid,'%s</%s>\n',levelTab,struName);
            else
                %insert template of rest of options' file
                fRid=fopen(fn.xml_template,'r');
                fseek(fRid,0,-1);
                while ~feof(fRid)
                    tline=fgetl(fRid);
                    fprintf(meta.fWid,'%s\n',tline);
                end
                fclose(fRid);
                fclose(meta.fWid);
            end
        end
        %when dne goin through, decrease level in 1
        level=level-1;
        
    end %of function [meta]= strprint( stru,meta)
    
end

function resampling(mouse, sess)

    fprintf('resampling\n')
    fprintf('============================================================\n')
    fn  = file_names(mouse, sess);
    q   = load(fn.ss_sess_info);
    rec = q.info.rec;
    tic;

    if ~exist(fn.fold_pr_sess, 'dir')
        mkdir(fn.fold_pr_sess)
    end
    
    for kr = 1:numel(rec)
        fprintf('============================================================\n')
        fprintf('rec : %s \n', rec(kr).name)

        fn = file_names(mouse, sess, rec(kr).name);
        fid   = fopen(fn.ss_rec, 'r');
        
        fs    = rec(kr).sampling_freq;
        fr    = 1000;
        bin   = floor(fs/2000);
        nChan = rec(kr).nChan;
                
        for ic = 1:nChan
            
            fprintf('chan: %3d  %7s  -> reading', ic, rec(kr).chan(ic).name)
            fseek(fid, (ic-1)*2, 'bof');
            
            [Y, n] = fread(fid, inf, 'int16', (nChan-1)*2);
            
            fprintf('   -> resampling')
            
            nb = floor(n/bin);
            nr = floor(nb/(fs/bin)*fr);
            
            Yb = mean(reshape(Y(1:bin*nb), bin, nb),1);
            
            Yr = interp1((1:nb)/(fs/bin), Yb, (1:nr)/fr);
            
            rsm.(rec(kr).chan(ic).name) = int16(Yr);
            fprintf('    %5d\n', round(toc))
        end
            
        fprintf('saving')
        save(fn.rsm_data, '-struct', 'rsm')
        fprintf('        done    %5d\n', round(toc))
        
        
        for ir = 1:numel(rec(kr).run)
            rec(kr).run(ir).start    = round(rec(kr).run(ir).offset/fs*fr);
            rec(kr).run(ir).duration = round(rec(kr).run(ir).nSampl/fs*fr);
        end
    end
    info = q.info;
    info.rec = rec;
    save(fn.sess_info, 'info')
     
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

function trial_prep(mouse, sess)

    fprintf('trial preparation\n')
    fprintf('============================================================\n')
    fn  = file_names(mouse, sess);
    q   = load(fn.sess_info);
    info = q.info;
    rec  = info.rec;
    
    for irec = 1:numel(rec)
        fn = file_names(mouse, sess, rec(irec).name);
        
        trial = struct();
        kt    = 0;
        
        for irun = 1:numel(rec(irec).run)
            fnam = fullfile(fn.fold_rd_sess, rec(irec).run(irun).behav_data);
            
            V = basicVoyerRead(fnam, 'sniff_ttl', 'sniff');
            trial = trial_build_01(trial);
            disp([irun, kt])
        end
        save(fn.sess_info, 'info')
        save(fn.trial,     'trial')
    end
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function tr = trial_build(tr)
        n_tr   = numel(V.trial);
        events = fields(V.event);
        start  = rec(irec).run(irun).start;

        for it = 1:n_tr
            if (it>1)&&(V.trial(it).trialstart < V.trial(it-1).trialstart)
                continue
            end

            kt = kt + 1;
            t  = V.trial(it);
            t0 = t.trialstart;
            tr(kt).run       = irun;
            tr(kt).tr_num    = it;
            tr(kt).start     = t0 + start;
            tr(kt).duration  = t.trialend - t0;
            tr(kt).odorTime  = t.fvOnTime*[1;1] + [0; t.fvdur] - t0*[1;1];
            tr(kt).odorName  = t.odor;
            tr(kt).odorConc  = t.odorconc;
            tr(kt).laserTime = t.laserontime*[1;1] + [0; t.duration_1] - t0*[1;1];
            tr(kt).laserAmpl = t.amplitude_1;
            tr(kt).stimID    = t.stimid;

            tr(kt).VoyerParameters = t;

            if it==n_tr
                t1 = V.event.breaks(2,end);
            else
                t1 = V.trial(it+1).fvOnTime-1;
            end

%             tr(kt).sniffWaveform = V.stream.sniff(t0:t1);

            for ie = 1:numel(events)
                in = (V.event.(events{ie})(2,:) > t0)&(V.event.(events{ie})(1,:)<t1);
                tr(kt).(events{ie}) = V.event.(events{ie})(in) - t0;
            end

        end
    end    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function tr = trial_build_01(tr)
        n_tr   = numel(V.trial);
        events = fields(V.event);
        start  = rec(irec).run(irun).start;
        for it = 1:n_tr
            if (it>1)&&(V.trial(it).trialstart < V.trial(it-1).trialstart)
                continue
            end
            kt = kt + 1;
            t  = V.trial(it);
            t0 = t.fvOnTime;
            tr(kt).run       = irun;
            tr(kt).tr_num    = it;
            tr(kt).start     = t0 + start;
            tr(kt).duration  = t.trialend - t0;
            tr(kt).odorTime  = t.fvOnTime*[1;1] + [0; t.fvdur] - t0*[1;1];
            tr(kt).odorName  = t.odor;
            tr(kt).odorConc  = t.odorconc;
            tr(kt).laserTime = t.laserontime*[1;1] + [0; t.duration_1] - t0*[1;1];
            tr(kt).laserAmpl = t.amplitude_1;
            tr(kt).stimID    = t.stimid;

            tr(kt).VoyerParameters = t;

            if it==n_tr
                t1 = V.event.breaks(2,end);
            else
                t1 = V.trial(it+1).fvOnTime-1;
            end
%             tr(kt).sniffWaveform = V.stream.sniff(t0:t1);

            for ie = 1:numel(events)
                in = (V.event.(events{ie})(2,:) > t0)&(V.event.(events{ie})(1,:)<t1);
                tr(kt).(events{ie}) = V.event.(events{ie})(in) - t0;
            end
        end
    end    
end

function V = basicVoyerRead(fnam, event_names, stream_names)

    if ~iscell(event_names)
        event_names = {event_names};
    end
    if ~iscell(stream_names)
        stream_names = {stream_names};
    end
    
    info = h5info(fnam);
    tr      = read_table(fnam);
    V.trial = tr;
    
  	events = read_events(fnam, event_names);
    for ke = 1:numel(event_names)
        V.event.(event_names{ke}) = events{ke};
    end
    
    [streams, breaks] = read_streams(fnam, stream_names);
    for ks = 1:numel(stream_names)
        V.stream.(stream_names{ks}) = streams(ks,:);
    end
    V.event.breaks = breaks;
    


    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function tr = read_table(fnam)
        
        n_tr = numel(info.Groups);
        table_name = info.Datasets.Name;
        table = h5read(fnam,['/', table_name]);
        ff = fieldnames(table);
        tr = struct();
        
        for kf = 1:numel(ff)
            for it = 1:n_tr
                if ischar(table.(ff{kf}))
                    tr(it).(ff{kf}) = deblank(table.(ff{kf})(:,it)');
                else
                    tr(it).(ff{kf}) = double(table.(ff{kf})(it,1));
                end
            end
        end        
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function [streams, breaks] = read_streams(fnam, sname)
        % read continious sniff waveform and check the packets interuption
        % streams{k} - continious waveforms. if there is a missing packet sn is
        % padded by zeros
        % breaks{k} - [2,nb] - sequence of miisng chnnks of data onset, offset
                
        disp('reading all streams...')
        n_str     = numel(sname);
        sample_on = zeros(1,1e7);
        streams   = zeros(n_str,1e7);
                
        for kt = 1:numel(info.Groups)
            group_name   = info.Groups(kt).Name;
            Events       = h5read(fnam, [group_name,'/Events']);
            
            for ks = 1:n_str
                streamPackets = h5read(fnam, [group_name,'/',sname{ks}]);
            
                for ie = 1:numel(Events.packet_sent_time)
                    nt  = int32(numel(streamPackets{ie}));
                    ind = Events.packet_sent_time(ie) - nt-1 +(1:nt);
                    if ind == 0
                        continue
                    end
                    streams(ks,ind) = streamPackets{ie};
                    sample_on(ind) = 1;
                end
            end
        end
        
        ind_end = ind(end)+10000;
        streams = streams(:,1:ind_end);
        
        breaks_on  = find(diff(sample_on)==-1);
        breaks_off = find(diff(sample_on)==1);
        
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
       
        breaks = [breaks_on(:)'; breaks_off(:)'];

    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function events = read_events(fnam, ename)
        % reads xx.h5 file trial-by-trail, extract continious sequence ofevents.
        % if the first event ==0, remove this event
        % create a matrix events [2, ne], where
        % ne - number of events
        % event(1,:), event(2,:) are onset and offess of the event
        
        disp('reading all events...')
        
        % read the initial trial's licks, see if there is an even or odd
        % number of licks. Licks are only transmitted when the beamstate is
        % off. If there is an odd number of lick times, then the first
        % value is off, if even then the first value is on. Want the first
        % value to be on, so discard the first value.
        

        n_names = numel(ename);
        events  = cell(1,n_names);
        for kn = 1:n_names
            events{kn} = zeros(1,1e5);
            ne        = 0;
       
            for kt = 1:numel(info.Groups)
                group_name  = info.Groups(kt).Name;
                eventPacket = h5read(fnam, [group_name,'/', ename{kn}]);
                event_tr    = cell2mat(eventPacket);
                ne_tr       = numel(event_tr);
                events{kn}(ne+(1:ne_tr)) = event_tr;
                ne = ne + ne_tr;
                if (kt == 1)&&mod(ne,2)
                    events{kn} = events{kn}(2:end);
                    ne    = ne-1;
                end
            end
            ne    = floor(ne/2);
            events{kn} = reshape(events{kn}(1:2*ne), 2, ne);
        end
    end
end

