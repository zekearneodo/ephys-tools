% =============  da_tools =================================================
% the set of tools for Data Analysis (da)
%                                           Dima Rinberg
%                                           version 1.01
%                                           date:   10.27.2009
%
%
% to activate this package call 
%   > da = da_tools();
% then use each individual function
% the list of function is below
% 
% fn = file_names(mouse, sess, disk) - the function creates a structure of file
%       and folder names correspondnet to this mouse and the session. disk
%       is the current working disk. Example:
%   > fn = da.file_names(553, 2, '/Volumes'RinbergLab')
%
%
% event_data_reformat(mouse, sess) - the function finds for a given mouse 
%       and session the program finds all txt even files and rewrite them 
%       in binary new format: 'uint32'. The old files are moved to archive
%       folder
%
%
% tr = read_behav_data(mouse, sess) - the function first cehck the presence of
%       behavioral '.mat' file for a given muse and session and reads it, 
%       if it does not find  '.mat' file it read '.txt' file and creates 
%       the strcuture 'tr'
%   > tr = da.read_behav_data(mouses, sess);
%
%
% event = read_event_data(mouse, sess, event_name) - the function reads 
%       binary event data for the mouse and session for the specific event 
%       file. The events are always goes in pairs, the first is 
%       trigger up, the second is down, the output of the function is an 
%       array (2,n_events). The events are time stamps in ms from the 
%       beginning of the recording
%   > event = da.read_event_data(mouse, sess, event_name);
%
%
% trig = read_trig_data(mouse, sess, trig_name) - the function reads 
%       binary event data for the mouse and session for the specific trigger 
%       file. Similar to events, but the trigger output is array
%       (1,n_trigger) of time stamps (in ms)
%   > trig = da.read_trig_data(mouse, sess, event_name);
%
%
% n_read = read_raw_data(mouses, sess, chan, offset, n_point) - the function 
%       reads the raw data file for agiven mouse, sess and chan,  offset is 
%       the number of points to skip from the beginning of the file, 
%       n_points is the number of point to read. The function returns n_read
%       the atual number of read points, and creates a global variable 
%       raw_data, where it store read raw data. This is done in order not 
%       to move large arrays of data. 
%   > global raw_data
%   > fs     = 32552
%   > n_read = da.read_raw_data(fnam, 100*fs, 10*fs);
%       the array raw data contains the data for 10 sec of recordings
%       starting 100 sec from the beginning of the file
%
%
% downsampling(mouse, sess, chan_name) - the function reads the raw data
%       from the chan file for a given mouse and session. (Chan_name is
%       specified as a string of only channel name: for example chan_name =
%       'Sniff', and the correspondent file fillbe 'rec01_Sniff.bin'.) The
%       function downsample it to 1 kHz sampling frequence and store it in
%       the file with the same name abut extension '.rsm'.
%   > da.downsampling(mouse, sess, 'Sniff')
%
%
% data = read_rsm_data(mouse, sess, chan_name) - read the
%       resampled data for mouse and sessinn for a given chan_name, returns
%       the data and number of points
%   > data = da.read_rsm_data(mouse, sess, 'Sniff');
%
% load_sess(mouse, sess) - the function creates a global structure d and
%       reads the following information for a given session: header.mat,
%       trais_mmm_ss.mat, spikes_mmm_ss.mat, and all resampled channels ,
%       like 'rec01_Sniff.rsm')
%   > d = load_sess(mouse, sess)
%
% build_event_file(mouse, sess) - the function reads the file 
%       'trials_mmm_ss.mat' and creates the file 'events_mmm_ss.mat', which
%       contains the strcuture 'event' with fields: laser_on, laser_off,
%       blank, odor_01, odor_02.... 
%
% all_mice() - create the list of all mice names which are availbel on the
%       curent  data disk
%
% all_sess(mouse) - create a list of all session numbers for a given mouse
%
% all_chan(mouse, sess) - creates a list of all channels, which are files with 
%       extension .bin for a given mouse and session
%
% [gf, gs] = subplot_fig(nx, ny, fig) -  creates a figure with subplots
%       nx - number of colomns, ny - number of raws, fig - the figure
%       handle, if fig does not exist, create a new figure. 
%       gf - figure handle, gs - an array [nx x ny] of subplot handles
%
% [rate, tt] = avearge_rate(sp, bin, limits) estimate psth from spike data
%       sp{1,n_tr} - a cell array of spike time for n_tr trials
%       bin        - bin size 
%       limits     - time axis limits, if no limits, min/max of sp
%
% [xx,yy] = raster(sp, limits) - builds arrays of x & y coordinates of spikes
%       for ratser plot
%       sp{1,n_tr} - cell of spike times, n_tr number of trials
%       limits     - time limits, if not specified defined by min & max od sp{},
%
% [rr, tt] = rate_matrix(sp, limits) - create a matrix rr(nt, n_tr) of zeros and ones for
%       spike response. nt - number of points in one trial, n_tr number of
%       trials. tt - is an optional output - time axis

function fcn = da_tools()
% da_tools provide a set of function for data analysis
% file_names - create a structure with file names
% read_raw_data
% read_behav_data
% read_event_data
% init_sess
% downsampling
% reformat even data

% commonly used global variables
% x - for large raw data arrays

fcn.file_names          = @file_names;
fcn.event_data_reformat = @event_data_reformat;
fcn.read_behav_data     = @read_behav_data;
fcn.read_event_data     = @read_event_data;
fcn.read_trigger_data   = @read_trigger_data;
fcn.read_raw_data       = @read_raw_data;
fcn.downsampling        = @downsampling;
fcn.read_rsm_data       = @read_rsm_data;
fcn.load_sess           = @load_sess;
fcn.build_event_file    = @build_event_file;
fcn.all_mice            = @all_mice;
fcn.all_sess            = @all_sess;
fcn.all_chan            = @all_chan;
fcn.average_rate        = @average_rate;
fcn.raster              = @raster;
fcn.rate_matrix         = @rate_matrix;


% fcn.read_spike_data     = @read_spike_data;
% fcn.init_sess           = @init_sess;
% fcn.spike_data_reformat = @spike_data_reformat;
% fcn.std_estimation      = @std_estimation;      
% fcn.viewer              = @viewer;
fcn.subplot_fig         = @subplot_fig;

end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function fn = file_names(mouse, sess)
global fn_disk


if exist('mouse', 'var')&&strcmp(mouse, 'ask')
    fn_disk = setup_disk(); fn = struct();
    return
end

if isempty(fn_disk)
    fn_disk = setup_disk();
end

if ismac
    fn.disk      = fullfile('/Volumes', fn_disk);
else
    fn.disk      = fn_disk;
end

fn.data_fold     = fullfile(fn.disk, 'Experiment', 'Data');
fn.analysis_fold = fullfile(fn.disk, 'Experiment', 'analysis');

if exist('mouse', 'var')
    if isfloat(mouse)
        mouse = sprintf('%03d', mouse);
    end
    fn.mouse_fold    = fullfile(fn.data_fold, sprintf('mouse_%s', mouse));
    fn.db_main_fold  = fullfile(fn.disk, 'Experiment', 'mouse_db');
end
    

if exist('sess', 'var')
    if ischar(sess)
        sess = str2double(sess);
    end

    fn.sess_fold     = fullfile(fn.mouse_fold, sprintf('session_%03d',sess));
    fn.ephys_data    = fullfile(fn.sess_fold, 'ephys data');
    fn.header        = fullfile(fn.ephys_data, 'header.mat');
    fn.behav         = fullfile(fn.sess_fold, sprintf('beh_%03d.txt', sess));
    fn.behav_mat     = fullfile(fn.sess_fold, sprintf('beh_%03d.mat', sess));
    fn.spikes        = fullfile(fn.analysis_fold, sprintf('spikes_%s_%02d.mat', mouse, sess));
    fn.trials        = fullfile(fn.analysis_fold, sprintf('trials_%s_%02d.mat', mouse, sess));
    fn.events        = fullfile(fn.analysis_fold, sprintf('events_%s_%02d.mat', mouse, sess));

    fn.sniff         = fullfile(fn.ephys_data, 'rec01_sniff.rsm');
    
    fn.db_sess_fold  = fullfile(fn.db_main_fold, sprintf('m%s_%02d', mouse, sess));
    fn.db_trials     = fullfile(fn.db_sess_fold, sprintf('trials_%s_%02d.mat', mouse, sess));
    fn.db_spikes     = fullfile(fn.db_sess_fold, sprintf('spikes_%s_%02d.mat', mouse, sess));
    fn.db_rsm_data   = fullfile(fn.db_sess_fold, sprintf('rsm_data_%s_%02d.mat', mouse, sess));
    fn.db_info       = fullfile(fn.db_sess_fold, sprintf('sess_info_%s_%02d.mat', mouse, sess));
    fn.db_sess       = fullfile(fn.db_sess_fold, sprintf('sess_%s_%02d.mat', mouse, sess));
    fn.db_events     = fullfile(fn.db_sess_fold, sprintf('events_%s_%02d.mat', mouse, sess));
    
    fn.sniff_waveforms = fullfile(fn.db_sess_fold, sprintf('sniff_%s_%02d.mat', mouse, sess)); 

end

    function disk = setup_disk()
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
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function event_data_reformat(mouse, sess_num)
% reformat event data files
% for a given mouse and session the program finds all txt even files and
% rewrite them in binary new format. The old files are moved to archive
% folder

    fn = file_names(mouse, sess_num);
    
    fold = fullfile(fn.ephys_data, 'archive');
    if ~exist(fold, 'dir')
        mkdir(fold)
    end
    
    
    list = dir(fn.ephys_data);

    for il = 1:length(list)
        fnam = list(il).name;
        if (~list(il).isdir)&&strcmp(fnam(end-4:end), '1.txt')
%             fnam_old = fullfile(fold, fnam);
%             if exist(fnam_old, 'file')
%                 delete(fnam_old)
%             end
            
            fnam_txt = fullfile(fn.ephys_data, fnam);
            disp(fnam_txt)
            a = textread(fnam_txt);
%             movefile(fnam_txt, fold);
        
            fnam_evn = fullfile(fn.ephys_data, sprintf('rec01_%s.evn',fnam(1:end-5)));
            disp(fnam_evn)
    
            fid = fopen(fnam_evn, 'w');
            fwrite(fid, a*1000, 'uint64');
            fclose(fid);
        end
    end

end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function tr = read_behav_data(mouse, sess)
% read the behavioral data for mouse and sess
% first the program checks if the correpondent '.mat' file exist, if yes,
% load the data from '.mat' file, if not, read the data from asqii file

% if the line starts with #p:, read it as a parameter and keep this
% parameter until it changes
% read the variable list. each variable line starts with '#v:'
% then read data into the data structure according to varibale list
% # - comments
% emty lines are ignored

fn = file_names(mouse, sess);

if exist(fn.behav_mat, 'file')
    disp(fn.behav_mat)
    q = load(fn.behav_mat);
    tr = q.tr;
    return
end

fid = fopen(fn.behav, 'r');
if fid < 0,
    warndlg(sprintf('file does not exist \n %s', fnam), 'Warning')
    return
end

% first read the variable list
not_eof = 1;    nv  = 0;
while not_eof
    tline       = fgetl(fid);
    not_eof = ~(isnumeric(tline)&&(tline == -1));
        
    if not_eof && strncmpi(tline, '#v:', 2)
        % variable line
        nv = nv + 1;
        [token, tline] = strtok(tline(3:end), ':');
        variable{nv} = strtrim(token);
        form{nv}     = strtrim(strtok(tline, ':'));
    end
end


% read data and parameters
fseek(fid, 0, 'bof');   
kt      = 0;   
tr      = [];     
not_eof = 1;
par     = struct;

while not_eof
    tline       = fgetl(fid);
    not_eof = ~(isnumeric(tline)&&(tline == -1));
    
    if not_eof && strncmpi(tline, '#p:',2)
        % parameter line
        [name, rem] = strtok(tline(3:end), ':');
        name        = strtrim(name);
        value       = strtok(rem(2:end), ':');
        if isnan(str2double(value))
            value = strtrim(value);
        else
            value = str2double(value);
        end
        
        % adding field to the structure
        par.(name) = value;
        param      = fields(par);
    end

    if not_eof&&(~isempty(tline))&&(tline(1) ~= '#')
        % data line
        kt  = kt + 1;
        for ip = 1:length(param)
            tr(kt).(param{ip}) = par.(param{ip});
        end
        
        for kv = 1:nv
            [val, tline] = strtok(tline);
            if isempty(findstr(form{kv}, 's'))
                tr(kt).(variable{kv}) = str2num(val);   % numeric
            else
                tr(kt).(variable{kv}) = strtrim(val);   % string
            end
        end                                                    
    end 
end

    
fclose(fid);

end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function ev = read_event_data(mouse, sess, event_name)
% read event data
% event data is stored as a sequence of 32-bit unsighed integer.
% each number is in ms from the beginning of the recording
    fn = file_names(mouse, sess);
    fnam = fullfile(fn.ephys_data, sprintf('rec01_%s.evn', event_name));
    if ~exist(fnam, 'file')
%         warndlg('file does not exist')
        ev = zeros(2,0);
        return
    end
    
    fid = fopen(fnam, 'r');
    [ev, ne]  = fread(fid, 'uint64');
    fclose(fid);
    if round(ne/2) == ne/2
        ev  = reshape(ev, 2, ne/2)/1000;
    else
        ev = reshape(ev(1:ne-1), 2, (ne-1)/2)/1000;
        disp('wrong number of events')
    end
    

end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function ev = read_trigger_data(mouse, sess, trig_name)
% read trigger data
% the same type of data as evn but the outcome consiste of single time
% stamps in ms
    fn = file_names(mouse, sess);
    fnam = fullfile(fn.ephys_data, sprintf('rec01_%s.evn', trig_name));
    if ~exist(fnam, 'file')
        warndlg('file does not exist')
        return
    end

    fid = fopen(fnam, 'r');
    ev  = fread(fid, 'uint64');
    ev = ev/1000;
    fclose(fid);
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function n_read = read_raw_data(mouses, sess, chan, offset, n_points)
% read data from the raw data file
% offset   - skip the first points from the beginning of the file
% n_points - number of points to read
% OUTPUT
% the result of the read function is stored in global variable raw_data
% n_read - actual number of read points

global raw_data

    fn = file_names(mouses, sess);
    fnam = fullfile(fn.ephys_data, sprintf('rec01_%s.bin', chan));
    if ~exist(fnam, 'file')
        warndlg('file does not exist')
        return
    end
    fid = fopen(fnam, 'r');
    fseek(fid, offset*2, 'bof');
    [raw_data, n_read] = fread(fid, n_points, 'int16');
    fclose(fid);

end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function downsampling(mouse, sess, chan_name)
% this function reads initial raw data and downsample it to 1 kHz sampling
% rate
% inputs: mouse, sess, rec, and chan_name
% the resampled data is stored in the file with the same name as raw data
% but with extension 'rsm'


global raw_data

tic
fs_rsm = 1000;             % downsampling frequency
bin    = 20;               % initial bin size for resizieng the raw data array
wind   = bin*1e6;         % the window size for initial data reading


fn = file_names(mouse, sess);
q  = load(fn.header);

fs     = q.header.rec.sampling;     % raw data sampling frequency
fs_bin = fs/bin;                    % sampling frequenct after rebinnging

nt_total = q.header.rec.number_of_samples+512;  % total number of points
nt_bin   = ceil(nt_total/bin);              % number of points in rebinned array

% reading raw data and rebining it
% fnam = fullfile(fn.ephys_data, sprintf('rec01_%s.bin', chan_name));
% disp(fnam)
% if ~exist(fnam, 'file')
%     warndlg('file does not exist')
%     return
% end

nt     = wind;
xb     = zeros(1, nt_bin);       % re-binned data
kt     = 0;                      % the beginning of raw data segment
kt_bin = 0;                      % the beginning of re-biibed data segment


while nt == wind
    nt  = read_raw_data(mouse, sess, chan_name, kt, wind);
    nt_bin = floor(nt/bin);
    nt     = nt_bin*bin;
    
    xb(kt_bin + [1:nt_bin]) = mean(reshape(raw_data(1:nt), bin, nt_bin),1);
    
    kt     = kt  + nt;
    kt_bin = kt_bin + nt_bin;

end

nt_bin = kt_bin;        % final numbe rof points in rebinned data array
xb = xb(1:nt_bin);

% resampling
nt_rsm = floor(nt_bin*fs_rsm/fs_bin);
xr     = interp1([1:nt_bin]/fs_bin, xb, [1:nt_rsm]/fs_rsm);

fnam = fullfile(fn.ephys_data, sprintf('rec01_%s.rsm', chan_name));

fid = fopen(fnam, 'w');
fwrite(fid, xr, 'int16');
fclose(fid);

toc

end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function s = read_rsm_data(mouse, sess, chan_name)
% read data from the rsm data file
% INPUT
% fnam - name of the file

    fn = file_names(mouse, sess);
    fnam = fullfile(fn.ephys_data, sprintf('rec01_%s.rsm', chan_name))
    if ~exist(fnam, 'file')
        warndlg('file does not exist')
        return
    end

    fid = fopen(fnam, 'r');
    [s, n_read] = fread(fid, 'int16');
    fclose(fid);
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function d = load_sess(mouse, sess_num)

    d.fn  = file_names(mouse, sess_num);
    q = load(d.fn.header);
    d.header = q.header;
    
    if exist(d.fn.trials, 'file')
        q = load(d.fn.trials);
        d.tr = q.tr;
    else
        d.tr = [];
    end
    
    if exist(d.fn.spikes, 'file')
        q = load(d.fn.spikes);
        d.unit = q.unit;
    else
        d.unit = [];
    end
    
    if exist(d.fn.sniff_waveforms, 'file')
        q = load(d.fn.sniff_waveforms);
        d.sniff = q.sniff;
    end
    
    list = dir(d.fn.ephys_data);

    for k = 1:length(list)
        fnam = list(k).name;
        if ~list(k).isdir
            if strcmp(fnam(end-3:end), '.rsm')
                [rn, s_name] = strtok(fnam, '_');
                s_name       = strtok(s_name(2:end), '.');
                d.(s_name)   = read_rsm_data(mouse, sess_num, s_name);
            end
            if strcmp(fnam(end-3:end), '.evn')
                [rn, e_name] = strtok(fnam, '_');
                e_name       = strtok(e_name(2:end), '.');
                d.(e_name)   = read_event_data(mouse, sess_num, e_name);
            end
        end
        
        
    end
        
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function build_event_file(mouse, sess)
    fn = file_names(mouse, sess);
    if ~exist(fn.trials, 'file')
        warndlg('trial file does not exist')
        return
    end
    q = load(fn.trials);
    tr = q.tr;
    event = struct();

    for it = 1:length(tr);
        if tr(it).laser_dur > 0
            add_event('laser_on',  tr(it).t.start + tr(it).t.laser(1));
            add_event('laser_off', tr(it).t.start + tr(it).t.laser(2));
        end
    
        if tr(it).final_valve_dur >0
            if tr(it).valve == 0
                name = 'blank';
            else
                name = sprintf('odor_%02d', tr(it).valve);
            end
            add_event(name, tr(it).t.start + tr(it).t.fv(1));
        end
    end

    event = orderfields(event);
    save(fn.events, 'event');

    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function add_event(name, t)
        if isfield(event, name)
            event.(name) = [event.(name), t];
        else
            event.(name)   = t;
        end
    end

end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function mm = all_mice()

    fn = file_names();
    list = dir(fn.data_fold); 
    nl = length(list);
    mm = cell(nl); km = 0;
    for il = 1:nl
        fnam = list(il).name;
        if ~isempty(findstr(fnam, 'mouse'))
            km = km + 1;
            mm{km} = fnam(7:end);
        end
    end
    mm = mm(1:km);
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function ss = all_sess(mouse)
    fn = file_names(mouse);
    list = dir(fn.mouse_fold); nl = length(list);
    ss = zeros(1,nl); ks = 0;
    for il = 1:nl
        fnam = list(il).name;
        if ~isempty(findstr(fnam, 'session'))
            ks = ks + 1;
            ss(ks) = str2double(fnam(9:end));
        end
    end
    ss = ss(1:ks);
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function chan = all_chan(mouse, sess)
    fn = file_names(mouse, sess);
    kc = 0;
    chan = cell(1,100);
    if exist(fn.ephys_data, 'dir');
        list = dir(fn.ephys_data);
        nl = length(list);
        for il = 1:nl
            fnam = list(il).name;
            if ~isempty(findstr(fnam, 'bin'))
                kc = kc + 1;
                chan{kc} = fnam(7:end-4);
            end
        end
    end
    chan = chan(1:kc);
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [rate, tt] = average_rate(sp, bin, limits)
    n_tr = length(sp);
    
    if exist('limits', 'var')
        t_mi = min(limits);
        t_ma = max(limits);
    else
        emp  = cellfun(@isempty, sp);
        t_mi = bin*floor(min(cellfun(@min, sp(~emp)))/bin)-bin;
        t_ma = bin* ceil(max(cellfun(@max, sp(~emp)))/bin)+bin;
    end
    tt   = t_mi+1:t_ma;
    nt   = length(tt);
    rate = zeros(1,nt);
    
    for it = 1:n_tr
        in = (sp{it} >  t_mi)&(sp{it} <= t_ma);
        sp_tr = ceil(sp{it}(in) - t_mi);
        rate(sp_tr) = rate(sp_tr) + 1;
    end
    
    nt1 = floor(nt/bin);
    rate = mean(reshape(rate(1:nt1*bin), bin, nt1),1)/n_tr*1000;
    tt   = mean(reshape(tt(1:nt1*bin),   bin, nt1),1);

end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [xx, yy] = raster(sp, varargin)

    
    n_arg = length(varargin);
    if n_arg/2 ~=round(n_arg/2)
        warning('wrong number or parameters !')
        return
    end
        
    v = struct;
    for iv = 1:2:n_arg
        v.(varargin{iv}) = varargin{iv+1};
    end
        
    if isfield(v , 'limits')
        t_mi = min(v.limits);
        t_ma = max(v.limits);
    else
        emp  = cellfun(@isempty, sp);
        t_mi = min(cellfun(@min, sp(~emp)));
        t_ma = max(cellfun(@max, sp(~emp)));
    end
    
    to_plot = 0;    if isfield(v, 'plot');    to_plot = v.plot;    end
    col     = 'k';  if isfield(v, 'color');   col =     v.color;   end
    gg      = gca;  if isfield(v, 'subplot'); gg  =     v.subplot; end
    

    n_tr = length(sp);
    
    if n_tr == 0
        warning('no spikes')
        return
    end
    
    xx = zeros(1,1e4);
    yy = zeros(1,1e4);
    n_sp = 0;
    for it = 1:n_tr
        nn = length(sp{it});
        xx(n_sp+(1:nn)) = sp{it};
        yy(n_sp+(1:nn)) = it*ones(1,nn);
        n_sp = n_sp + nn;
    end
    
    in_lim = (xx(1:n_sp)>t_mi)&(xx(1:n_sp) < t_ma);
    xx = xx(in_lim);
    yy = yy(in_lim);

%     subplot(gg) 
%     plot(xx, yy, '.', 'MarkerSize', 8, 'Color', col)
%     axis([t_mi, t_ma, 0, it+1]);
%     
%     if nargout > 0
%         varargout{1} = xx;
%         varargout{2} = yy;
%     end

end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [rr, varargout] = rate_matrix(sp, limits)
    if exist('limits', 'var')
        t_mi = min(limits);
        t_ma = max(limits);
    else
        emp  = cellfun(@isempty, sp);
        t_mi = min(cellfun(@min, sp(~emp)));
        t_ma = max(cellfun(@max, sp(~emp)));
    end

    n_tr = length(sp);
    n_tp = t_ma - t_mi + 1;
    rr   = zeros(n_tp, n_tr);
    for jt = 1:n_tr
        in = (sp{jt} > t_mi)&(sp{jt} <= t_ma);
        rr(ceil(sp{jt}(in)-t_mi(1)),jt) = 1;
    end
    
    if nargout > 0
        varargout = t_mi:t_ma;
    end
end








% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function sp = read_spike_data(fnam)
% read event data
% event data is stored as a sequence of 32-bit unsighed integer.
% each number is in ms from the beginning of the recording
fid = fopen(fnam, 'r');
sp  = fread(fid, 'uint32');
fclose(fid);
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function init_sess(mouse, snum)

fn = file_names(mouse, snum);
disp(fn)

load(fn.header)

if ~exist(fn.analysis_fold, 'dir')
    mkdir(fn.analysis_fold)
end

sess.mouse = mouse;
sess.num   = snum;
sess.date  = header.date;
sess.time  = header.time;

for ir = 1:length(header.rec)
    sess.rec(ir).start_time        = header.rec(ir).start_time;
    sess.rec(ir).sampling          = header.rec(ir).sampling;
    sess.rec(ir).number_of_samples = header.rec(ir).number_of_samples;
    for ic = 1:length(header.rec(ir).Achan)
        
                              name    = header.rec(ir).Achan(ic).name;
        sess.rec(ir).chan(ic).name    = name;
        sess.rec(ir).chan(ic).gain    = header.rec(ir).Achan(ic).preamp_gain;
        sess.rec(ir).chan(ic).range   = header.rec(ir).Achan(ic).input_range;
        sess.rec(ir).chan(ic).lp_filt = header.rec(ir).Achan(ic).high_cut_freq;
        sess.rec(ir).chan(ic).hp_filt = header.rec(ir).Achan(ic).low_cut_freq;
        
        if strncmpi(name, 'csc', 3)
            fnam = sprintf('rec%02d_S%02d.bin', ir, str2double(name(4:end)));
        else
            fnam = sprintf('rec%02d_%s.bin', ir, name);
        end
        sess.rec(ir).chan(ic).rd_fnam = fnam;
    end
    
    for ic = 1:length(header.rec(ir).Dchan);
                               name = header.rec(ir).Dchan(ic).name;
        sess.rec(ir).event(ic).name = name;
        sess.rec(ir).event(ic).fnam  = sprintf('rec%02d_%s.evn', ir, name);
    end
end
    
save(fn.sess_info, 'sess');

end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function spike_data_reformat(mouse, sess_num, unit)
    fn = file_names(mouse, sess_num);
    load(fn.sess_info);
    start_t = sess.rec(1).start_time;

    [fnam, path] = uigetfile([fn.sess_fold, '/*.txt']);
    fnam = fullfile(path, fnam);
    fprintf(fnam)
    fprintf('\n')

    sp = textread(fnam);
%     sp = round((sp-start_t)/1000);
    sp = round((sp-sp(1))/1000);

    fnam = fullfile(fn.analysis_fold, sprintf('rec01_unit%02d.spk', unit));
    fid = fopen(fnam, 'w');
    fwrite(fid, sp, 'uint32');
    fclose(fid);
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function std_estimation(mouse, sess_num, chan)
    global x

    filt = [300, 5000];

    fn = fcn.file_names(mouse, sess_num);
    load(fn.sess_info)
    fs = sess.rec(1).sampling/1000;

    wind   = 10e3;
    nchan  = length(chan);
    m_chan = max(chan);
    [bf, af] = butter(4, filt/(fs*1000/2));

    kt0 = 0;
    np_read = round(wind*fs);
    np      = np_read;

    Y        = zeros(np, m_chan);
    np_total = 0;
    Z1       = zeros(1, m_chan);
    Z2       = zeros(1, m_chan);

    while np == np_read,
        for ic = chan
            fnam = fullfile(fn.ephys_data, sess.rec(1).chan(ic).rd_fnam);
            np   = read_raw_data(fnam, kt0, np_read);
            Y([1:np],ic) = x;
        end
        Y = Y([1:np],:);
        Y = filtfilt(bf, af, Y);
        Y = Y - mean(Y,2)*ones(1, size(Y,2));
        Y = Y - ones(np,1)*mean(Y,1);
        np_total = np_total + np; 
        Z2       = Z2 + sum(Y.^2, 1);
        kt0      = kt0 + np;
    
        fprintf(sprintf('%d  %d \n', round(kt0/(fs*1000)), mean(Z2/np_total)))
    end

    Ystd = sqrt(Z2/np_total);

    for ic = chan,
        sess. rec(1).chan(ic).Ystd = Ystd(ic);
    end

    save(fn.sess_info, 'sess');
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function viewer(t0, wind, mouse, sess_num, chan, bckg, filt, unit, event)
global x
% input arguments:
% t0 - initial time [ms]
% wind - window     [ ms]
% mouse
% sess_num 
% chan - the list of chan, if chan is the cell array, plot as individual terodes
% bckg - background channel ('all')
% filt - filter [Hz]
% unit - spike unit numbers, if cell array plot units for each tetrode
% event - list of events ('all')

if iscell(chan)
    tetr = chan;
    chan = cell2mat(chan);
else
    for ic = 1:length(chan)
        tetr{ic} = chan(ic);
    end
end

if ischar(bckg)&&strcmpi(bckg, 'all')
    bckg = chan;
end

disp(chan)
disp(bckg)

all_chan = unique([chan, bckg]);
nchan    = length(chan);
m_chan   = max(all_chan);

fn = file_names(mouse, sess_num);
% load(fn.sess_info)
load(fn.header)
r = header.rec;
fs = r.sampling;
% fs = r.number_of_samples/(r.last_timeStamp - r.start_time)*1000;
% fs = 3.255208332432239e+04/1000;

kt0    = floor(t0*fs);
nt     = floor(wind*fs);

Y      = zeros(nt, m_chan);

% read raw data
for ic = all_chan
    fnam = fullfile(fn.ephys_data, sprintf('rec01_S%02d.bin', ic));
%     fnam = fullfile(fn.ephys_data, [header.rec(1).Achan(ic).binary_name, '.bin'])
    disp(fnam)
    nt   = read_raw_data(fnam, kt0, nt);
    Y(:,ic) = x;
end
    
% filtering
if exist('filt', 'var')&&(~isempty(filt))
    [bf, af] = butter(4, filt/(fs/2));
    Y(:,all_chan) = filtfilt(bf, af, Y(:,all_chan));
end

% background subtruction
if ~isempty(bckg)
    Ybckg = mean(Y(:,bckg), 2);
    Y(:,chan) = Y(:,chan) - Ybckg*ones(1,nchan);
end

% offset estimation
Ystd   = mean(std(Y(:,chan)));
offset = 20*Ystd;

% plotting
figure(4), clf
col    = 'bgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmyk';
% col    = 'bkmbkmbkmbkmbkmbkmbkmbkmbkmbkmbkmbkmbkmbkm';
% col    = 'kbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbk';

gs = subplot(1,1,1);
set(gs, 'position', [0.05, 0.05, 0.9, 0.9], ...
    'XLim',         t0+[0,wind], ...
    'YTick',        [1:nchan], ....
    'YGrid',        'on', ...
    'XGrid',         'on', ...
    'NextPlot',     'add')

kc = 0;
for it = 1:length(tetr);
    y0{it} = [];
    for ic = tetr{it};
        kc = kc + 1;
        plot(t0+[1:nt]/fs, kc + Y(:,ic)/offset, col(it)),
        YTick{kc} = num2str(ic);
        y0{it} = [y0{it}, kc]; 
    end
end

% spike plotting
if exist('unit', 'var')&&(~isempty(unit))
    
    for iu = 1:length(unit);
        fnam = fullfile(fn.analysis_fold, sprintf('spikes_%s.spk', unit{iu}));
        disp(fnam)
        sp = read_spike_data(fnam);
        sp = sp/1000;
        in = find((sp > t0)&(sp < t0 + wind));
        if ~isempty(in)
            plot(ones(nchan,1)*sp(in)', (y0{it}-0.3-iu*0.05)'*ones(1,length(in)), ...
                ['x',col(iu)], 'MarkerSize', 8, ...
                'LineWidth', 2);    %,   'Color',      col{iu})
        end
        
    end
end

% event plotting
if exist('event', 'var')&&(~isempty(event))
    for ie = 1:length(event),
        fnam = fullfile(fn.ephys_data, sprintf('rec01_%s.evn', event{ie}));
%         fnam = fullfile(fn.analysis_fold, sess.rec(1).event(ie).fnam);
        disp(fnam)
        ev = read_event_data(fnam);
        ev = ev/1000;

        in = find((ev(2,:) > t0)&(ev(1,:) < (t0+wind)));
        int = ev(:,in);
        n_ev = length(in);
        if n_ev > 0,
            int(1,1)    = max([int(1,1), t0]);
            int(2,n_ev) = min([int(2,n_ev), t0+wind]);
        end

        plot(int, (-ie/5)*ones(size(int)), 'LineWidth', 4, 'Color', col(ie)), hold on
    end
end

set(gs, 'Ylim', [-1,nchan+1], ....
    'YTickLabel',   YTick, ...
    'NextPlot',     'replace');



end
