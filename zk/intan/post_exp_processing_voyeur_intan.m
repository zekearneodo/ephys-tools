% Data preparation for spike sorting
% Dima Rinberg, August 2013
%
% ZK - Mar 2014:
% Based on post_exp_processing_voyeur_ndm. (in current)
% if prog_ephys is Intan, it reads channels in the intan way.
% As of March 5, it is ready and tested for intan with data on channels
% from:
%   -amplifier
%   -board_adc
%   -board_dig_in
% Will need update if other channels (with other sampling freqs) will be
% used.
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

function post_exp_processing_voyeur_intan(mouse, sess,stat)
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
pp.ch_correct     = @chan_info_correct;

if ~exist('mouse', 'var')    
    mouse = 'ZKanesthM71';
    sess  = '069';
    %stat = '01';
end

if nargin>2 && ~isempty(stat) && ~strcmp(stat,'local')
    pull_data(mouse,sess,stat);
end

read_info(mouse, sess);
ss_prep(mouse, sess)
%xml_make(mouse,sess);
resampling(mouse, sess)
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
fprintf('\n *** Making Info File ***\n\n');

fn = file_names(mouse, sess)

[info, rec] = read_log_file();


if strcmp(info.progr_ephys,'SpikeGL')
    fprintf('Recording sytem is SpikeGL\n')
    ftype  = {'ephys_data', 'ephys_meta', 'behav_data'};
    rec   = read_data_dir(rec);
    rec   = read_meta_data_SGL(rec);
else
    fprintf('Recording sytem is Intan\n')
    ftype  = {'ephys_data', 'behav_data'};
    rec   = read_data_dir(rec);
    rec   = read_meta_data_intan(rec);
end
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
fprintf('\nMetadata and sess info read\n\n');

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
        fprintf('Reading data dir...')
        % the program reads raw data sess folder, sort files by records and runs,
        % and reads meta information for records and runs
        % the program uses the 'rec' structure saved in 'ss_data\se_xxx_xx' folder
        % it uses the paprmaters from these structure
        % if it finds a record with a new name it uses the paarmeters of the last
        % record.
        % the output is saved back to the same file, same structure
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
                    %warnglg('new record name')
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
        fprintf(' done.\n')
    end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rec = read_meta_data_SGL(rec)
        fprintf('--- reading meta data for spikeGL recording system \n');
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
            in_digital = find(strcmp({chan_rd.type}, 'D')&~strcmpi({chan_rd.name}, 'pltrig'));
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
            iPltrig=find(strcmpi({chan_rd.name}, 'pltrig'));
            rec(kr).pltrig_num_rd  = chan_rd(iPltrig).num;
            
            offset = 0;
            for ir = 1:numel(rec(kr).run),  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                fn_meta = fullfile(fn.fold_rd_sess, rec(kr).run(ir).ephys_meta);
                meta   = read_meta_data(fn_meta);
                nSampl = meta.fileSizeBytes/meta.nChans/2;
                
                rec(kr).run(ir).nSampl = nSampl;
                rec(kr).run(ir).offset = offset;
                offset = offset + nSampl;
                
            end,    % ir ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        end,    % kr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
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
    end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rec = read_meta_data_intan(rec)
        %read the first segment of each rec to obtain the
        %frequency_parameters from the Intan file.
        fprintf('*Reading meta data...\n')
        chanAccount={'amplifier',...
            'board_adc',...
            'aux_input',...
            'board_dig_in',...
            'board_dig_out',...
            'supply_voltage',...
            'temp_sensor'};
        chanCodes={'e',...
            'c',...
            'a',...
            'd',...
            'o',...
            'v',...
            't'};
        keepFields={'notes',...
            'frequency_parameters'};
        %channels I know how to handle (for now), that will go to the .dat
        %file
        datChans={'amplifier',...
            'board_adc',...
            'board_dig_in',...
            };
        datFactors=[5 32768/5 5];
        
        for kr = 1:numel(rec)
            fprintf('Rec %s\n',rec(kr).name)
            fprintf('-\n')
            fnData   = fullfile(fn.fold_rd_sess, rec(kr).run(1).ephys_data);
            dataFiles = dir(fullfile(fnData,'*.rhd'));
            firstFile = fullfile(fnData,dataFiles(1).name);
            fprintf('Reading first rhd file to get metadata...',rec(kr).name)
            intMeta   = read_intan_chunk(firstFile,1);
            fprintf(' done\n',rec(kr).name)
            intanMeta = struct; %structure to keep track of intan complete meta (w.o. the big data)
            
            %total number of channels (total, total)
            %count the files
            nChan_rd=0;
            for chanType=chanAccount;
                nChan_rd=nChan_rd+intMeta.(['num_' chanType{1} '_channels']);
                if intMeta.(['num_' chanType{1} '_channels'])>0
                    intanMeta.([chanType{1} '_channels'])=intMeta.([chanType{1} '_channels']);
                end
            end
            
            for inField=keepFields
                intanMeta.(inField{1})=intMeta.(inField{1});
            end
            
            rec(kr).nChan_rd      = nChan_rd;
            rec(kr).sampling_freq = intMeta.frequency_parameters.amplifier_sample_rate;
            %the gain values are fixed for the Intan RHD200
            rec(kr).gain          = 200;
            rec(kr).rangeMax      = 2500;
            rec(kr).rangeMin      = -2500;
            rec(kr).intanMeta     = intanMeta;
            
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
                chCode=char(chan_rd(ic).num);
                type = chCode(1);
                iCh  = str2num(chCode(2:end));
                
                chan(kc).num_rd = chCode;
                chan(kc).type_rd=char(chanAccount(strcmp(type,chanCodes)));
                chan(kc).tidx_rd=iCh;
                chan(kc).factor =datFactors(strcmp(chan(kc).type_rd,datChans));
                
                chan(kc).name   = chan_rd(ic).name;
                if isfield(chan_rd,'group')
                    chan(kc).group  = chan_rd(ic).group;
                end
            end
            
            %remove the channels that I don't still know how to handle (or
            %tha will go to a different file)
            dataChans=zeros(1,numel(chan));
            for iDT=1:numel(datChans);
                dataChans=dataChans | ismember({chan.type_rd},datChans{iDT});
            end
            rec(kr).nChan = sum(dataChans);
            rec(kr).chan  = chan(dataChans);
            
            %warn that there were channels I don't know what to do
            if rec(kr).nChan<nChan
                warning('There were channels of a type I dont know what to do about in rec %s',rec(kr).name)
            end
            % digital channels
            in_digital = find(strcmp({chan_rd.type}, 'D')&~strcmpi({chan_rd.name}, 'pltrig'));
            nDigChan   = numel(in_digital);
            if nDigChan > 0
                for kc = 1:nDigChan
                    ic = in_digital(kc);
                    digChan(kc).name   = chan_rd(ic).name;
                    
                    chCode=chan_rd(ic).num;
                    type = chCode(1);
                    iCh  = str2num(chCode(2:end));
                
                    digChan(kc).num_rd = chCode;
                    digChan(kc).type_rd=chanAccount(strcmp(type,chanCodes));
                    digCchan(kc).tidx_rd=iCh;
                end
                rec(kr).digChan  = digChan;
            end
            rec(kr).nDigChan = nDigChan;
            
            % PL trigger
            iPltrig=find(strcmpi({chan_rd.name}, 'PLtrig'));
            rec(kr).pltrig_num_rd  = chan_rd(iPltrig).num;
            
            for krun=1:numel(rec(kr).run)
                rec(kr).run(krun).nSampl = 0;
                rec(kr).run(krun).offset = 0;
            end
        end
    end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function ch = read_chan_config(fnam)
        disp(['reading channel config file ' fnam '...']);
        if ~exist(fnam, 'file')
            error('no channel_config file')
        end
        fid  = fopen(fnam, 'r');
        if strcmpi(info.progr_ephys,'SpikeGL')
            str = textscan(fid, '%d  %s   %s   %s');
            num = double(str{1});
        elseif strcmpi(info.progr_ephys,'Intan')
            str = textscan(fid, '%s  %s   %s   %s');
            num = str{1};
        else
            warning('Unknown ephys recoding program');
            str = textscan(fid, '%d  %s   %s   %s');
            num = double(str{1});
        end
        type = str{2};
        name = str{3};
        numel(name);
        group= str{4};
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
  
function chan_info_correct(mouse, sess)
% the program read the following inofmration:
% 1 - log file
% 2 - raw data file directory
% 3 - electrode channel config
% 4 - meta data for individual recordings
% It prints the structure rec for debugging.
% It saves the data into ss_fold, sess_info.mat file
global info;
fprintf('\n *** Making Info File ***\n\n');

fn = file_names(mouse, sess)

[info, rec] = read_log_file();


if strcmp(info.progr_ephys,'SpikeGL')
    fprintf('Recording sytem is SpikeGL\n')
    ftype  = {'ephys_data', 'ephys_meta', 'behav_data'};
    rec   = read_data_dir(rec);
    rec   = read_meta_data_SGL(rec);
else
    fprintf('Recording sytem is Intan\n')
    ftype  = {'ephys_data', 'behav_data'};
    rec   = read_data_dir(rec);
    rec   = read_meta_data_intan(rec);
end
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
fprintf('\nMetadata and sess info read\n\n');

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
        fprintf('Reading data dir...')
        % the program reads raw data sess folder, sort files by records and runs,
        % and reads meta information for records and runs
        % the program uses the 'rec' structure saved in 'ss_data\se_xxx_xx' folder
        % it uses the paprmaters from these structure
        % if it finds a record with a new name it uses the paarmeters of the last
        % record.
        % the output is saved back to the same file, same structure
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
                    %warnglg('new record name')
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
        fprintf(' done.\n')
    end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rec = read_meta_data_SGL(rec)
        fprintf('--- reading meta data for spikeGL recording system \n');
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
            in_digital = find(strcmp({chan_rd.type}, 'D')&~strcmpi({chan_rd.name}, 'pltrig'));
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
            iPltrig=find(strcmpi({chan_rd.name}, 'pltrig'));
            rec(kr).pltrig_num_rd  = chan_rd(iPltrig).num;
            
            offset = 0;
            for ir = 1:numel(rec(kr).run),  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                fn_meta = fullfile(fn.fold_rd_sess, rec(kr).run(ir).ephys_meta);
                meta   = read_meta_data(fn_meta);
                nSampl = meta.fileSizeBytes/meta.nChans/2;
                
                rec(kr).run(ir).nSampl = nSampl;
                rec(kr).run(ir).offset = offset;
                offset = offset + nSampl;
                
            end,    % ir ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        end,    % kr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
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
    end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rec = read_meta_data_intan(rec)

        for kr = 1:numel(rec)
            fprintf('Rec %s\n',rec(kr).name)
            fprintf('-\n')
            fnData   = fullfile(fn.fold_rd_sess, rec(kr).run(1).ephys_data);
            dataFiles = dir(fullfile(fnData,'*.rhd'));
            firstFile = fullfile(fnData,dataFiles(1).name);
            fprintf('Reading first rhd file to get metadata...',rec(kr).name)
            intMeta   = read_intan_chunk(firstFile,1);
            fprintf(' done\n',rec(kr).name)
            intanMeta = struct; %structure to keep track of intan complete meta (w.o. the big data)
            
            %total number of channels (total, total)
            %count the files
            nChan_rd=0;
            for chanType=chanAccount;
                nChan_rd=nChan_rd+intMeta.(['num_' chanType{1} '_channels']);
                if intMeta.(['num_' chanType{1} '_channels'])>0
                    intanMeta.([chanType{1} '_channels'])=intMeta.([chanType{1} '_channels']);
                end
            end
            
            for inField=keepFields
                intanMeta.(inField{1})=intMeta.(inField{1});
            end
            
            rec(kr).nChan_rd      = nChan_rd;
            rec(kr).sampling_freq = intMeta.frequency_parameters.amplifier_sample_rate;
            %the gain values are fixed for the Intan RHD200
            rec(kr).gain          = 200;
            rec(kr).rangeMax      = 2500;
            rec(kr).rangeMin      = -2500;
            rec(kr).intanMeta     = intanMeta;
            
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
                chCode=char(chan_rd(ic).num);
                type = chCode(1);
                iCh  = str2num(chCode(2:end));
                
                chan(kc).num_rd = chCode;
                chan(kc).type_rd=char(chanAccount(strcmp(type,chanCodes)));
                chan(kc).tidx_rd=iCh;
                chan(kc).factor =datFactors(strcmp(chan(kc).type_rd,datChans));
                
                chan(kc).name   = chan_rd(ic).name;
                if isfield(chan_rd,'group')
                    chan(kc).group  = chan_rd(ic).group;
                end
            end
            
            %remove the channels that I don't still know how to handle (or
            %tha will go to a different file)
            dataChans=zeros(1,numel(chan));
            for iDT=1:numel(datChans);
                dataChans=dataChans | ismember({chan.type_rd},datChans{iDT});
            end
            rec(kr).nChan = sum(dataChans);
            rec(kr).chan  = chan(dataChans);
            
            %warn that there were channels I don't know what to do
            if rec(kr).nChan<nChan
                warning('There were channels of a type I dont know what to do about in rec %s',rec(kr).name)
            end
            % digital channels
            in_digital = find(strcmp({chan_rd.type}, 'D')&~strcmpi({chan_rd.name}, 'pltrig'));
            nDigChan   = numel(in_digital);
            if nDigChan > 0
                for kc = 1:nDigChan
                    ic = in_digital(kc);
                    digChan(kc).name   = chan_rd(ic).name;
                    
                    chCode=chan_rd(ic).num;
                    type = chCode(1);
                    iCh  = str2num(chCode(2:end));
                
                    digChan(kc).num_rd = chCode;
                    digChan(kc).type_rd=chanAccount(strcmp(type,chanCodes));
                    digCchan(kc).tidx_rd=iCh;
                end
                rec(kr).digChan  = digChan;
            end
            rec(kr).nDigChan = nDigChan;
            
            % PL trigger
            iPltrig=find(strcmpi({chan_rd.name}, 'PLtrig'));
            rec(kr).pltrig_num_rd  = chan_rd(iPltrig).num;
            
        end
    end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function ch = read_chan_config(fnam)
        disp(['reading channel config file ' fnam '...']);
        if ~exist(fnam, 'file')
            error('no channel_config file')
        end
        fid  = fopen(fnam, 'r');
        if strcmpi(info.progr_ephys,'SpikeGL')
            str = textscan(fid, '%d  %s   %s   %s');
            num = double(str{1});
        elseif strcmpi(info.progr_ephys,'Intan')
            str = textscan(fid, '%s  %s   %s   %s');
            num = str{1};
        else
            warning('Unknown ephys recoding program');
            str = textscan(fid, '%d  %s   %s   %s');
            num = double(str{1});
        end
        type = str{2};
        name = str{3};
        numel(name);
        group= str{4};
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

function ss_prep(mouse, sess)
    
    fn = file_names(mouse, sess);
    if ~exist(fn.ss_sess_info, 'file')
        error('no sess_info file')
    end
    
    if isnumeric(mouse)
        mouse=num2str(mouse);
    end
    
    fprintf('*************************************************\n')
    fprintf(' PRE-PROCESSING SESSION  %s OF MOUSE %s \n', sess, mouse)
    fprintf('**************************************************\n')
    tic;    
    q = load(fn.ss_sess_info);
    rec = q.info.rec;
    progr_ephys=q.info.progr_ephys;

    for kr = 1:numel(rec),   % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fprintf('==================================================\n')
        fprintf('rec: %s \n', rec(kr).name)
        
        ipl      = rec(kr).pltrig_num_rd;         % PL trigegr channel number(or code)
        chan     = rec(kr).chan;                  % channel parameters
        nChan    = rec(kr).nChan;                 % number of channels
        nChan_rd = rec(kr).nChan_rd;              % number of raw data channels
        
        pl_extr = strcmp(rec(kr).pl_trigger, 'yes')&~isempty(ipl);
        
        % create a record SS folder and SS file
        fn = file_names(mouse, sess, rec(kr).name);
        mkdir(fn.fold_ss_rec)
        fid_rec = fopen(fn.ss_rec, 'w');
        fseek(fid_rec,0,'bof');
        runOffset=0;
        
        for ir = 1:numel(rec(kr).run)
            fprintf('--------------------------------------------------\n')
            fprintf('run: %d ', rec(kr).run(ir).num)
            
            run    = rec(kr).run(ir);                  % run parameters
            
            % open raw data file
            fn_data = fullfile(fn.fold_rd_sess, run.ephys_data);
            
            % do the run getting depending on the recording system
            fun=sprintf('%s_get_run_data', lower(progr_ephys));
            getRunDataFcn=eval(['@' fun ';']);
            
           
            %this function gets the channels, removes the pl trigger
            %and
            %saves the channel where it corresponds in the final binary
            %file.
            feval(getRunDataFcn,fn_data);
            %to save back the run parameters, if any, that are defined or
            %modified in getRunDataFcn
            rec(kr).run(ir)=run;

        end
        fclose(fid_rec);
    
    end
    
    q.info.rec=rec;
    save(fn.ss_sess_info,'-struct','q','info');
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function spikegl_get_run_data(fn_data)
        %the whole data for the run is in one single file
        fid_r = fopen(fn_data, 'r');
        
        fprintf('file: %s\n ', fn_data)
        %get the pl_trigger
        if pl_extr
            fprintf('pl_trigger \n')
            pl      = read_pl_trigger(fid_r, ipl, nChan_rd);
            pl_extr = pl_check(pl)
        end
        
        %get all the channels and remove the pl_trigger
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
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    function intan_get_run_data(fn_data)
        %the data for the run is in several one-minute chunks, in the
        %folder.
        %get the list of files
        chunkOffset=0;
        fprintf(' (folder: %s)\n', fn_data)
        fl=dir(fullfile(fn_data,'*.rhd'));
        [~,iflSort]=sort({fl.name});
        fl=fl(iflSort);
        run.offset=runOffset;  
        run.nSampl=0;
        for ifl=1:numel(fl)
            fprintf(' Chunk %3d/%d: %s ',ifl,numel(fl), fl(ifl).name)
            %get the chunk
            %
            fprintf('read.. ')
            chunk=read_intan_chunk(fullfile(fn_data,fl(ifl).name),1);
            %reserve the memory for the matrix that will then be written
            chunkDatSize=length(chunk.amplifier_data);
            dataMatrix=zeros(chunkDatSize,nChan);
            
            %get the pl_trigger
            if pl_extr
                fprintf('getting plTrig.. ')
                pl      = extract_intan_chan(ipl, chunk);
                pl_extr = pl_check(pl)
            end
            
            %get all the other channels
            fprintf('get %4d ch.. ', nChan)
            for kc = 1:nChan,
                X = extract_intan_chan(kc, chunk);
                X = X - mean(X);
                
                if pl_extr
                    X = remove_pl_trigger(X, pl);
                end
                
                dataMatrix(:,kc)=int16(round(X*chan(kc).factor));
            end
            chunkOffset=chunkOffset+length(dataMatrix);
            fprintf('write.. ')
            fwrite(fid_rec,transpose(dataMatrix),'int16');
            run.nSampl = run.nSampl+length(dataMatrix);
            fprintf('(%d s)\n',round(toc));
        end
        runOffset  = runOffset+run.nSampl;
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
    function pl_extr = pl_check(pl)
        if isempty(pl)
            pl_extr = 0;
        else
            pl_freq = rec(kr).sampling_freq/mean(diff(pl));
            if (pl_freq < 50)||(pl_freq > 70)
                pl_extr = 0;
            end
        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function Y = read_analog_channel(fid, ich, nch)
        fseek(fid, 2*(ich-1), 'bof');
        Y = fread(fid, inf, 'int16', (nch-1)*2);
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function Y=extract_intan_chan(kc, chunk)
        %extract the channel corresponding to chCode from the intan
        %structure.
        %chCode is the num_rd.
        
        %find out the channel index (in the channel structure).
        
        type = char(chan(kc).type_rd);
        iCh  = chan(kc).tidx_rd;
        
        %for now assume there were no lost samples; go ahead and read.
        Y=chunk.([type '_data'])(iCh,:);
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
        return
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
            fprintf('    (%5d s)\n', round(toc))
        end
            
        fprintf('saving')
        save(fn.rsm_data, '-struct', 'rsm')
        fprintf('        done    (%5d s)\n', round(toc))
        
        
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
