% data management tools:
%   Dima Rinberg
%   Ezequiel Arneodo
%                       last updaate:   Nov 18, 2013
%
% collection of the programs for data management in the lab
% version 031 is the same as 03, but instead of using the 'Laser' signal,
% uses the 'TrPin'.
%
% ---------------------------------------------------------------------------------------
% file_names - creates the structure fn for all file names 
% fn = file_names()
% fn = file_names(mouse, sess)
% fn = file_names(mouse_sess,rec);
%
% ---------------------------------------------------------------------------------------
% read_raw_data_info(mouse, sess)
% the program read the following inofmration:
% 1 - log file
% 2 - raw data file directory
% 3 - electrode channel config
% 4 - meta data for individual recordings
% It prints the structure rec for debugging.
% It saves the data into ss_fold, sess_info.mat file
%
% ---------------------------------------------------------------------------------------
% sniffs_to_trial(mouse,sess)
%
% order of post_spike functions:
%  -sniff_analysis_04_zk
%  -afta_da_sorting(outside, have to include here)
%  -trial_prep
%  -trial_shift_estim
%  -trial_alignment 
%  -sniffs_to_trials
%  -spikes_to_trials

% mouse = 'KPawakeM72';
% sess = '001';
% rec = 'a';

function dm = data_management_tools_032_kp()
global dm

    dm.file_names         = @file_names;
    dm.read_raw_data_info = @read_raw_data_info;
    dm.trial_prep         = @trial_prep;
    dm.basicVoyerRead     = @basicVoyerRead;
    dm.trial_shift_estim  = @trial_shift_estim;
    dm.trial_alignment    = @trial_alignment;
    dm.sniffs_to_trials   = @sniffs_to_trials;
    dm.spikes_to_trials   = @spikes_to_trials;
    
    dm.subplot_fig        = @subplot_fig;
    
% pp.read_data_dir  = @read_data_dir;
% pp.ss_prep        = @ss_prep;
% pp.resampling     = @resampling;
end

function fn = file_names(mouse, sess, rec, stat)
%global fn_disk

%     if exist('mouse', 'var')&&strcmp(mouse, 'ask')
%         fn_disk = setup_disk(); 
%         fn = struct();
%     end
% 
%     if isempty(fn_disk)
%         fn_disk = setup_disk();
%     end

    local_disk  = '/Volumes/spikefolder';
%     local_disk = '/Users/kpenikis/Documents/Rinberg/Experiment';
%     local_disk  = '/home/zeke/data';
%     stat_dist   = '/stations';
    
    fn.disk = local_disk;
    fn.fold_rd_data  = fullfile(fn.disk, 'raw_data');
    fn.fold_ss_data  = fullfile(fn.disk, 'ss_data');
    fn.fold_pr_data  = fullfile(fn.disk, 'pr_data');
    fn.fold_config   = fullfile(fn.disk, 'SpikeGL_config');
    
%     fn.fold_sd_data  = fullfile(stat_disk, '');
    if nargin > 1
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
        
        fn.fold_ss_rec  = fullfile(fn.fold_ss_sess,  sprintf('rec_%s',rec));
        fn.ss_rec       = fullfile(fn.fold_ss_rec,   sprintf('rec_%s.dat',rec));
        fn.rsm_data     = fullfile(fn.fold_pr_sess,  sprintf('%s_%s_%s_rsm.mat', mouse, sess, rec));
        fn.trials       = fullfile(fn.fold_pr_sess,  sprintf('%s_%s_%s_trial.mat', mouse, sess, rec));
        fn.spikes       = fullfile(fn.fold_pr_sess,  sprintf('%s_%s_%s_spikes.mat', mouse, sess, rec));
        fn.sniffs       = fullfile(fn.fold_pr_sess,  sprintf('%s_%s_%s_sniff.mat', mouse, sess, rec));
    end
    
    if nargin >3
        if isnumeric(stat)
            stat = sprintf('%02d', stat);
        end
        
        fn.fold_sd_data  = fullfile(stat_disk, sprintf('stat_%s',stat),'raw_data');
        fn.fold_sd_mouse = fullfile(fn.fold_sd_data,   sprintf('mouse_%s',mouse));
        fn.fold_sd_sess  = fullfile(fn.fold_sd_mouse,  sprintf('sess_%s',sess));
    end
    

    

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function disk = setup_disk()
        
        disks = {'/Volumes/Exp', '/experiment', '/Volumes/spikefolder','/home/zeke/data'};
        
        [sel, ok] = listdlg('Liststring', disks, ...
            'SelectionMode',    'single', ...
            'ListSize',         [150, 50], ...
            'PromptString', 'Choose disk:');
        disk = disks{sel}; 
        fprintf('disk: %s \n', disk)
    end

    function disk = setup_disk1()
        if ismac
            list = dir('/Volumes');
            kd = 0; d = cell(0);
            for il = 1:length(list)
                if list(il).name(1) ~='.'
                    kd = kd + 1;
                    d{kd}= list(il).name;
                end
            end
        else
            import java.io.*; 
            f=File('');
            r=f.listRoots;
            for i=1:numel(r)
                d{i} = char(r(i));
            end
        end
        disp(d)
 
        [sel, ok] = listdlg('Liststring', d, ...
            'SelectionMode',    'single', ...
            'ListSize',         [150, 50], ...
            'PromptString', 'Choose disk:');
        disk = d{sel}; 
        fprintf('disk: %s \n', disk)
    end

end

function read_raw_data_info(mouse, sess)
% the program read the following inofmration:
% 1 - log file
% 2 - raw data file directory
% 3 - electrode channel config
% 4 - meta data for individual recordings
% It prints the structure rec for debugging.
% It saves the data into ss_fold, sess_info.mat file

    fn = file_names(mouse, sess);

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
        save('recdbg.mat','rec')
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
                chan(kc).name    = chan_rd(ic).name;
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
            rec(kr).pltrig_num_rd  = find(strcmp({chan_rd.name}, 'PLtrig'));

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
            if ~exist(fnam, 'file')
                error('no channel_config file')
            end

            fid  = fopen(fnam, 'r');
            str = textscan(fid, '%d  %s   %s');
            num = double(str{1});
            type = str{2};
            name = str{3};
            fclose(fid);

            ch = struct();
            for kc = 1:numel(num)
                ch(kc).num  = num(kc);
                ch(kc).name = name{kc};
                ch(kc).type = type{kc};
            end
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
% reads the voyeur parameters
%builds the trial structure with that data
global trial
    fprintf('trial preparation\n')
    fprintf('============================================================\n')
% basically brings up the info struct and gives nickname to rec info sub-array
    fn  = file_names(mouse, sess);
    q   = load(fn.sess_info);
    info = q.info;
    rec  = info.rec;
    
    for irec = 1:numel(rec)
        fn = file_names(mouse, sess, rec(irec).name); % pulls up all files with that mouse,sess,rec in the name
        
        trial = struct();
        kt    = 0;
        
        for irun = 1:numel(rec(irec).run)
            fnam = fullfile(fn.fold_rd_sess, rec(irec).run(irun).behav_data);
            fprintf('Rec %s, run %d\n',rec(irec).name,irun);
            V = basicVoyerRead(fnam, 'sniff_ttl', 'sniff');
            trial = trial_build_04(trial);
            disp([irun, kt])
        end
        save(fn.sess_info, 'info')
        save(fn.trials,     'trial')
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

            tr(kt).VoyeurParameters = t;

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
%         events = fields(V.event);
        runStart  = rec(irec).run(irun).start;
        v_events  = {'sniff_ttl', 'breaks'};
        t_events  = {'sniffTTLTimes', 'breakTimes'};
        laser_on = [1 1 1 1 0 1 1 0 1];

        for it = 1:n_tr

            kt = kt + 1;
            t  = V.trial(it);
            t0 = t.fvOnTime;
            tr(kt).run           = irun;
            tr(kt).runTrialNum   = it;
            tr(kt).runTrialStart = t0;
            tr(kt).recTrialStart = t0+runStart;
            tr(kt).runTrialDur   = t.trialend - t0;
            
            tr(kt).start      = [];             % will be filled after alignment 
            tr(kt).duration   = [];
            
            tr(kt).odorTimes  = t.fvOnTime*[1;1] + [0; t.fvdur] - t0*[1;1];
            tr(kt).odorName   = t.odor;
            tr(kt).odorConc   = t.odorconc;
            
            tr(kt).laserTimes = t.laserontime*[1;1] + [0; laser_on(irun)*t.duration_1/1000] - t0*[1;1];
            tr(kt).laserAmp   = t.amplitude_1*laser_on(irun);
            tr(kt).laserPower = [];
            
            tr(kt).stimID     = t.stimid;

            tr(kt).VoyeurParameters = t;
            

            if t.trialdur>0
                t1 = t0 - 2000;
                t2 = t0 + 6000;

                for ie = 1:numel(v_events)
                    in = (V.event.(v_events{ie})(2,:) > t1)&(V.event.(v_events{ie})(1,:) < t2);
                    tr(kt).(t_events{ie}) = V.event.(v_events{ie})(:,in) - t0;
                end
            end

        end
    end

    function tr = trial_build_02(tr)
        %written for the october 2013 versions of behavior programs of
        %chronic1 rig.

        n_tr   = numel(V.trial);
        events = fields(V.event);
        runStart  = rec(irec).run(irun).start;
        v_events  = {'sniff_ttl', 'breaks'};
        t_events  = {'sniffTTLTimes', 'breakTimes'};
        laser_on = [1 0 1 1 0];
        %laser_on = [1 1 1 1 1];
        for it = 1:n_tr

            kt = kt + 1;
            t  = V.trial(it);
            t0 = t.trialstart;
            tr(kt).run           = irun;
            tr(kt).runTrialNum   = it;
            tr(kt).runTrialStart = t0;
            tr(kt).recTrialStart = t0+runStart;
            tr(kt).runTrialDur   = t.trialend - t0;
            
            tr(kt).start      = [];             % will be filled after alignment 
            tr(kt).duration   = [];
            
            tr(kt).odorTimes  = t.fvOnTime*[1;1] + [0; t.fvdur] - t0*[1;1];
            tr(kt).odorName   = t.odor;
            tr(kt).odorConc   = t.odorconc;
            
            tr(kt).laserTimes = t.laserontime*[1;1] + [0; laser_on(irun)*t.duration_1/1000] - t0*[1;1];
            tr(kt).laserAmp   = t.amplitude_1*laser_on(irun);
            tr(kt).laserPower = [];
            
            tr(kt).stimID     = t.stimid;

            tr(kt).VoyeurParameters = t;
            

            if t.trialdur>0
                t1 = t0 - 2000;
                t2 = t0 + 6000;

                for ie = 1:numel(v_events)
                    in = (V.event.(v_events{ie})(2,:) > t1)&(V.event.(v_events{ie})(1,:) < t2);
                    tr(kt).(t_events{ie}) = V.event.(v_events{ie})(:,in) - t0;
                end
            end
            
            %Reaconditionning of the trials to:
            %When no odor was presented, odorTimes on and off are equal, and
            %odorConc=0.
            if strcmp(t.odor,'dummy')
                tr(kt).odorTimes(2,:) = tr(kt).odorTimes(1,:);
                tr(kt).odorName   = t.odor;
                tr(kt).odorConc   = t.odorconc;
            end
            
            

        end
    end

    function tr = trial_build_03(tr)
        %written for the december 2013 versions of behavior programs of
        %chronic1 rig.

        n_tr   = numel(V.trial);
        events = fields(V.event);
        runStart  = rec(irec).run(irun).start;
        v_events  = {'sniff_ttl', 'breaks'};
        t_events  = {'sniffTTLTimes', 'breakTimes'};
        laser_on  = zeros(1,numel(rec(irec).run));
        for it = 1:n_tr

            kt = kt + 1;
            t  = V.trial(it);
            t0 = t.trialstart; % start trial from beginning of beh file (run)
            tr(kt).run           = irun; 
            tr(kt).runTrialNum   = it;
            tr(kt).runTrialStart = t0;
            tr(kt).recTrialStart = t0+runStart;
            tr(kt).runTrialDur   = t.trialend - t0;
            
            tr(kt).start      = [];             % will be filled after alignment 
            tr(kt).duration   = [];
            
            tr(kt).odorTimes  = t.fvOnTime*[1;1] + [0; t.fvdur] - t0*[1;1]; 
            tr(kt).odorName   = t.odor;
            tr(kt).odorConc   = t.odorconc;
            
            tr(kt).laserTimes = t.laserontime*[1;1] + [0; laser_on(irun)*t.duration_1/1000] - t0*[1;1];
            tr(kt).laserAmp   = t.amplitude_1*laser_on(irun);
            tr(kt).laserPower = [];
            
            tr(kt).stimID     = t.stimid;

            tr(kt).VoyeurParameters = t;
            

            if t.trialdur>0
                t1 = t0 - 2000;
                t2 = t0 + 6000;

                for ie = 1:numel(v_events)
                    in = (V.event.(v_events{ie})(2,:) > t1)&(V.event.(v_events{ie})(1,:) < t2);
                    tr(kt).(t_events{ie}) = V.event.(v_events{ie})(:,in) - t0;
                end
            end
            
            %Reaconditionning of the trials to:
            %When no odor was presented, odorTimes on and off are equal, and
            %odorConc=0.
            if strcmp(t.odor,'dummy')
                tr(kt).odorTimes(2,:) = tr(kt).odorTimes(1,:);
                tr(kt).odorName   = t.odor;
                tr(kt).odorConc   = t.odorconc;
%             %When there was odor, laser hast to be set off.
%             else
%                 laser_on(irun)    = 0;
%                 tr(kt).laserTimes = t.laserontime*[1;1] + [0; laser_on(irun)*t.duration_1/1000] - t0*[1;1];
%                 tr(kt).laserAmp   = t.amplitude_1*laser_on(irun);
            end
        end
    end

    function tr = trial_build_04(tr)
        %written for the jan 2014 protocols of acute1 rig.

        n_tr   = numel(V.trial);
        events = fields(V.event);
        runStart  = rec(irec).run(irun).start;
        v_events  = {'sniff_ttl', 'breaks'};
        t_events  = {'sniffTTLTimes', 'breakTimes'};
        laser_on  = zeros(1,numel(rec(irec).run));
        for it = 1:n_tr

            kt = kt + 1;
            t  = V.trial(it);
            t0 = t.trialstart;
            tr(kt).run           = irun;
            tr(kt).runTrialNum   = it;
            tr(kt).runTrialStart = t0;
            tr(kt).recTrialStart = t0+runStart; %trial start from the beginning of the rec
            tr(kt).runTrialDur   = t.trialend - t0;
            
            tr(kt).start      = [];             % will be filled after alignment 
            tr(kt).duration   = [];
            
            tr(kt).odorTimes  = t.fvOnTime*[1;1] + [0; t.fvdur] - t0*[1;1];
            tr(kt).odorName   = t.odor;
            tr(kt).odorConc   = t.odorconc;
            
            tr(kt).laserTimes = t.laserontime*[1;1] + [0; laser_on(irun)*t.duration_1/1000] - t0*[1;1];
            tr(kt).laserAmp   = t.amplitude_1*laser_on(irun);
            tr(kt).laserPower = [];
            
            tr(kt).stimID     = t.stimid;

            tr(kt).VoyeurParameters = t;
            

            if t.trialdur>0
                t1 = t0 - 2000;
                t2 = t0 + 6000;

                for ie = 1:numel(v_events)
                    in = (V.event.(v_events{ie})(2,:) > t1)&(V.event.(v_events{ie})(1,:) < t2);
                    tr(kt).(t_events{ie}) = V.event.(v_events{ie})(:,in) - t0;
                end
            end
            
            %Reaconditionning of the trials to:
            %When no odor was presented, odorTimes on and off are equal, and
            %odorConc=0.
            if strcmp(t.odor,'dummy')
                tr(kt).odorTimes(2,:) = tr(kt).odorTimes(1,:);
                tr(kt).odorName   = t.odor;
                tr(kt).odorConc   = t.odorconc;
%             %When there was odor, laser hast to be set off.
%             else
%                 laser_on(irun)    = 0;
%                 tr(kt).laserTimes = t.laserontime*[1;1] + [0; laser_on(irun)*t.duration_1/1000] - t0*[1;1];
%                 tr(kt).laserAmp   = t.amplitude_1*laser_on(irun);
            end
        end
    end

end

function V = basicVoyerRead(fnam, event_names, stream_names)
%reads the h5 files generated by voyeur.
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

function trial_shift_estim(mouse, sess, rec)
%estimates the misalingment (in trials) of the h5 file of timestapms
%respect to the ephys recordings
%It was made for trigger signals in the beh files that count the time since
%the beginning of the trial.
%For use with the Trial Pin and back compatibility,adapted that the time tt
%does not add any other delay if 'runTrialStart' and 'TrPin' are the
%signals.
%global info

    fn = file_names(mouse, sess);
    q = load(fn.sess_info);
    info = q.info;
    
    
    if ~exist('rec', 'var')
        reclist = {q.info.rec.name};
    elseif ~iscell('rec')
        reclist = {rec};
    end   
% recorded channels to compare trials from: ephys and behavior
    eTrig_chan = 'TrPin';
    bTrig_chan = 'runTrialStart';
    trial_shift = -3:3; nts = numel(trial_shift);
    
    for irec = 1:numel(reclist)
        nrec = find(strncmpi(cellstr({q.info.rec.name}),reclist(irec),1))
        nrun = numel(info.rec(nrec).run);
        fn = file_names(mouse, sess, reclist{irec});
        q = load(fn.trials);
        trial = q.trial;
        
        q = load(fn.rsm_data, eTrig_chan); %load Tr_pin channel only from rsm.mat file


        % ephys trial start trigger
        eSync = find(diff(q.(eTrig_chan) > 1e3)==1);  %start times of trials as recorded in TrPin by SpikeGL
        
        % beh trial start trig
        if ~ strcmp('runTrialStart',bTrig_chan)
            tt = [1;1]*[trial.runTrialStart] + [trial.(bTrig_chan)];
            bSync = tt(1,:);  %some default if old recording?
        else
            bSync=[trial.runTrialStart];  %start times of trials as recorded by Voyeur (see fct trial_prep)
                                                %NOTE: these times are wrt
                                                %beginning of record
        end
        
        figure(11), clf
        
        for irun = 1:nrun
            % e-phys trigger
            run_start = info.rec(nrec).run(irun).start;
            run_dur   = info.rec(nrec).run(irun).duration;
    
            in = (eSync >= run_start)&(eSync <= run_start + run_dur);
            eSyncRun = eSync(in);  %just the trials in current run
            
            % beh trigger
            in = [trial.run] == irun;
            bSyncRun = bSync(in);
            inGood   = [trial(in).runTrialDur]>0;
    
            % estimation synchronization error for each trial shift value
            error = zeros(1,nts);
            alpha = zeros(1,nts);
            err_min = 1e10;
            
            for it = 1:nts       % nts is the number of shifts to try
                % defining the shifted trigger times
                if trial_shift(it) >= 0    %for positive trial_shift values, eSyncRun is shifted
                    eTimes = eSyncRun((1+trial_shift(it)):end);
                    bTimes = bSyncRun;
                    ind    = inGood;    % beh trials with duration >0 ... why?
                else    %for negative trial_shift values, bSyncRun is shifted
                    eTimes = eSyncRun;
                    bTimes = bSyncRun((1-trial_shift(it)):end);
                    ind    = inGood((1-trial_shift(it)):end);
                end
                
                ne = numel(eTimes);
                nb = numel(bTimes);
                nn = min([ne, nb]);

                %  and allinging triggers by the first trial
                eTimes = eTimes(2:nn) - eTimes(1);
                bTimes = bTimes(2:nn) - bTimes(1);
                ind    = ind(2:nn);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                alpha(it) = mean(eTimes(ind)./bTimes(ind));    
                error(it) = sqrt(mean((eTimes(ind) - alpha(it)*bTimes(ind)).^2));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                subplot(nrun,3,3*irun)
                plot(trial_shift(it), error(it), 'o', 'Color', 'r', 'LineWidth', 2)
                hold on
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if error(it)<err_min    %check if this shift new min error
                    err_min = error(it);
                    imin    = it;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    subplot(nrun,3,3*irun-2)
                    plot(eTimes(ind), bTimes(ind), 'x')   
                    title(trial_shift(it))
         
                    subplot(nrun,3,3*irun-1)
                    plot(eTimes(ind) -alpha(it)*bTimes(ind), 'x')
                end
                drawnow
            end
            
            
            info.rec(nrec).run(irun).trialShift          = trial_shift(imin)
            info.rec(nrec).run(irun).frequencyCorrection = alpha(imin);
            disp(info.rec(nrec).run(irun))
        end
    end
    save(fn.sess_info, 'info')
    
end

function trial_alignment(mouse, sess, rec)
%aligns the trials and the ephys file.
    fn = file_names(mouse, sess);
    q = load(fn.sess_info);
    info = q.info;
    
    q.info.rec.name
%     indrec = find(q.info.rec.name=={rec})
    
    if ~exist('rec', 'var')
        rec = {q.info.rec.name};
    elseif ~iscell('rec')
        rec = {rec};
    end   
    


    for irec = 1:numel(rec)       % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fn = file_names(mouse, sess, rec{irec});
        
        eTrig_chan = 'TrPin';
        
        
        % loading ephys sync channel
        q  = load(fn.rsm_data, eTrig_chan);
        eSync = find(diff(q.(eTrig_chan) > 1e3)==1);
        spikegl_trials = length(eSync)
        
        % loading behavioral trials
        q = load(fn.trials);
        trial = q.trial;
        voyeur_trials = length(trial)
        irun  = 0;
        
        for it = 1:numel(trial),   % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % check if it is a new run, if it is, get run data for all sync triggers
          
            if trial(it).run ~= irun
                irun = trial(it).run;
                run  = info.rec(irec).run(irun)

                fprintf('rec %s, run %d: ===========================\n', info.rec(irec).name, irun)
                
                % defined the ephys sync triggers
                eSyncRun = eSync((eSync > run.start)&(eSync < run.start+run.duration)) - run.start;
                
                % make shift correction
                if isfield(run, 'trialShift')&&~isempty(run.trialShift)
                    trialShift = info.rec(irec).run(irun).trialShift;
                else
                    trialShift = 0;
                end                    
                   
                if isfield(run, 'frequencyCorrection')&&~isempty(run.frequencyCorrection)
                    freqCorr = run.frequencyCorrection;
                else
                     freqCorr = 1;
                end
                
                if trialShift >= 0
                    eFirst = 1 + trialShift;
                    bFirst = 1;
                else
                    eFirst = 1;
                    bFirst = 1 - trialShift;
                end
                
                % shift the beginning of the run to the first trial
                eStart   = run.start + eSyncRun(eFirst);
                eSyncRun = eSyncRun  - eSyncRun(eFirst);
                
                % find first behavioral trial for a given run
                ii = find(([trial.run]==irun)&([trial.runTrialNum] == bFirst), 1);
                bStart = trial(ii).runTrialStart;
                
                %find the average cumulative shift per trial
                bSync = [trial.runTrialStart];
                nn    = min(numel(eSyncRun),numel(trial));
                drift = mean(diff(eSync(1:nn))-diff(bSync(1:nn)));
                
            end
            
            t0 = trial(it).runTrialStart;
            
            % find the start of e-phys trial in the 1000 ms vicinity
            % around behavioral trial
            ii = find(abs(eSyncRun - freqCorr*(t0 - bStart))<3500, 1);
            if ~isempty(ii)
                trial(it).start = eSyncRun(ii) + eStart - trial(it).laserTimes(1)*strcmp(eTrig_chan,'Laser');
            else
                trial(it).start = 0;
            end
            
            fprintf('%d    %d \n', it, trial(it).start)
        end
        
%         save(fn.trials, 'trial')
        
    end
    
    


end

function sniffs_to_trials(mouse, sess, rec)
% the function reads the strucutre sniff.mat from xxx_sniff.mat file and
% sort the sniff parameters into individula trials, save the results tback
% to trial structure to the file xxx_trials.mat
% the mosue and sess must be specifiied. 
% rec is optional
% if rec is specified, it sort sniffs to trials for a specific recors,
% otehrwise, it reads the lsit of all records for a given session

global trial

    fn = file_names(mouse, sess);
    q = load(fn.sess_info);
    
    if ~exist('rec', 'var')
        rec = {q.info.rec.name};
    elseif ~iscell('rec')
        rec = {rec};
    end   

    for ir = 1:numel(rec)
        fn = file_names(mouse, sess, rec{ir});
        q = load(fn.trials);
        trial = q.trial;
        
        q = load(fn.sniffs);
        sniff = q.sniff;
        n_sn  = numel(sniff);
        
        sn_t_zer     = ones(3,1)*[sniff.t0] + [sniff.t_zer];
        sn_t_zer_fit = ones(4,1)*[sniff.t0] + reshape([sniff.t_zer_fit], 4, n_sn);
        
        for it = 1:numel(trial)
            if isempty(trial(it).start)||(trial(it).start == 0)
                continue
            end
            
            t0 = trial(it).start;
            t1 = t0 - 2000;
            t2 = t0 + 6000;
            
            in_sn = (sn_t_zer(3,:) > t1)&(sn_t_zer(1,:) < t2);
            
            trial(it).sniffZeroTimes      = sn_t_zer([1,2], in_sn) - t0;
            trial(it).sniffParabZeroTimes = sn_t_zer_fit(:,in_sn)  - t0;
            
            fprintf('%d    %d \n', it, sum(in_sn))
            
        end
        
        save(fn.trials, 'trial')
    end
    
    

end

function spikes_to_trials(mouse, sess, rec)

global trial

    fn = file_names(mouse, sess);
    q = load(fn.sess_info);
    
    if ~exist('rec', 'var')
        rec = {q.info.rec.name};
    elseif ~iscell('rec')
        rec = {rec};
    end   

    for ir = 1:numel(rec)
        fn = file_names(mouse, sess, rec{ir});
        q = load(fn.trials);
        trial = q.trial;
        
        q = load(fn.spikes);
        
        for iu = 1:numel(q.unit)
            sp{iu} = q.unit(iu).times;
        end
        
        
        for it = 1:numel(trial)
            if isempty(trial(it).start)||(trial(it).start == 0)
                continue
            end
            
            t0 = trial(it).start;
            t1 = t0 - 2000;
            t2 = t0 + 6000;
            
            fprintf('%d    ', it)
            for iu = 1:numel(q.unit)
                trial(it).spikeTimes{iu} = sp{iu}((sp{iu}>t1)&(sp{iu}<t2)) - t0;
                fprintf('%2i   ', numel(trial(it).spikeTimes{iu}));
            end
            
            fprintf('\n')
        end
        
    end
    
    save(fn.trials, 'trial');

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

