%this function makes the file structure for one element of ephys data.
%receives parameters of the recording, and returns a structure with file
%and folder names.
%Oct 2013, DR
%Modified by ZK

%It should be placed in /basefolder/ephysDataManagement/current/include
%functions calling it should include that folder in the path by calling
%addpath('/basefolder/ephysDataManagement/current/include') 
%or something like that.

function fn = file_names(mouse, sess, rec, stat)
%this function makes the file structure for one element of ephys data.
%receives parameters of the recording, and returns a structure with file
%and folder names.
persistent local_disk;

%if it is lunux or mac, checks the hostname
[~,computerName]=system('hostname');
%check if it is the computing server and do everything that is particular
%to it
if strcmpi(strtrim(computerName),'flipper')
        %the file server
        config_fold = '/usr/local/kluster/config';
        cf=fopen(fullfile(config_fold,'experiment_folder'),'r');
        local_disk=fgetl(cf);
        fclose(cf);
        
        fn.xml_template= fullfile(config_fold,'template.xml');
        fn.ndm_def_par = fullfile(local_disk,'ndm_defaults');
        stat_disk='/stations';
end

%look or prompt for local_disk
if isempty(local_disk)
%sets up base folder for the experiment
%if its the computational server, load the configurations.
%careful when you replace this, keep in mind it will affect all the
%computers.
%Here, of all places, changes MUST be backwards compatible
fprintf('Setting base folder for computer is %s)\n',computerName);
switch strtrim(computerName)
    case 'flipper'
        cf=fopen(fullfile(config_fold,'experiment_folder'),'r');
        local_disk=fgetl(cf);
        fclose(cf);
    case 'Dima-MacBookPro4.local'
        %zk notebook
        local_disk  = '/Users/zeke/experiment';
    case 'kristina-macpro.ersp.nyumc.org';
        local_disk = '/Volumes/spikefolder';
    case 'kristinikissmbp.vz30.nyumc.org';
%         local_disk = '/Users/kpenikis/Documents/Rinberg/Experiment';
        local_disk = '/Volumes/spikefolder';
    otherwise
        %assume any other with the server mounted
        local_disk = '/Volumes/spikefolder';
end %switch computerName

        
%check if the location exits, otherwise warn and prompt
if ~exist(local_disk,'dir')
    warning('Folder %s not found!! (is it mounted)?\n Prompting location\n',local_disk);
    local_disk = uigetdir('','Select base folder for your ephys');
end

fprintf('Base folder set to %s \n',local_disk);
end

%fprintf('Configuring file_names with base folder %s \n',local_disk);
fn.base_folder   = local_disk;
fn.fold_rd_data  = fullfile(local_disk, 'raw_data');
fn.fold_ss_data  = fullfile(local_disk, 'ss_data');
fn.fold_pr_data  = fullfile(local_disk, 'pr_data');
fn.fold_an_data  = fullfile(local_disk, 'an_data');
fn.fold_exp_data = fullfile(local_disk, 'export_data');
fn.fold_unit_db  = fullfile(local_disk, 'units');
fn.fold_st_meta  = fullfile(local_disk, 'stimuli');
fn.fold_config   = fullfile(local_disk, 'SpikeGL_config');
fn.fold_prb      = fullfile(local_disk, 'probe_definitions');

if exist('stat_disk','var')
    fn.fold_sd_data  = fullfile(stat_disk, 'var');
end


if nargin > 0
    if isnumeric(mouse)
        mouse = sprintf('%04d', mouse);
    end
    fn.mouse = mouse;
    
    if nargin>1
        if isnumeric(sess)
            sess = sprintf('%03d', sess);
        end

        fn.sess  = sess;
        fn.fold_rd_mouse = fullfile(fn.fold_rd_data,  sprintf('mouse_%s', mouse));
        fn.fold_rd_sess  = fullfile(fn.fold_rd_mouse, sprintf('sess_%s', sess));
        fn.log           = fullfile(fn.fold_rd_sess,  sprintf('log_%s_%s.txt', mouse, sess));
        
        fn.fold_ss_sess  = fullfile(fn.fold_ss_data,  sprintf('ss_%s_%s', mouse, sess));
        fn.fold_pr_sess  = fullfile(fn.fold_pr_data,  sprintf('m%s_%s', mouse, sess));
        fn.ss_sess_info  = fullfile(fn.fold_ss_sess,  sprintf('%s_%s_sess_info.mat', mouse, sess));
        fn.sess_info     = fullfile(fn.fold_pr_sess,  sprintf('%s_%s_sess_info.mat', mouse, sess));
        fn.clInfo_file   = fullfile(fn.fold_pr_sess,  sprintf('%s_%s_cl_info.mat',mouse,sess));
    end
end

%if rec was specified
if nargin >2 && ~isempty(rec)
    if isnumeric(rec)
        rec = sprintf('%02d', rec);
    end
    
    fn.fold_ss_rec  = fullfile(fn.fold_ss_sess,   sprintf('rec_%s',rec));
    fn.rd_rec_bn    = sprintf('%s',rec);
  
    fn.ss_rec       = fullfile(fn.fold_ss_rec,    sprintf('rec_%s.dat',rec));
    fn.ss_kk2_prm   = fullfile(fn.fold_ss_sess,   sprintf('rec_%s.prm', rec));
    fn.ss_xml       = fullfile(fn.fold_ss_sess,   sprintf('rec_%s.xml',rec));
    fn.lfp          = fullfile(fn.fold_ss_rec,    sprintf('rec_%s.lfp',rec));
    fn.rsm_data     = fullfile(fn.fold_pr_sess,   sprintf('%s_%s_%s_rsm.mat', mouse, sess, rec));
    fn.trial        = fullfile(fn.fold_pr_sess,   sprintf('%s_%s_%s_trial.mat', mouse, sess, rec));
    fn.spikes       = fullfile(fn.fold_pr_sess,   sprintf('%s_%s_%s_spikes.mat', mouse, sess, rec));
    fn.sniffs       = fullfile(fn.fold_pr_sess,   sprintf('%s_%s_%s_sniff.mat', mouse, sess, rec));
    fn.evt_m        = fullfile(fn.fold_pr_sess,   sprintf('%s_%s_%s_event.mat', mouse, sess, rec));
    fn.evt_h        = fullfile(fn.fold_pr_sess,   sprintf('%s_%s_%s_event.h5', mouse, sess, rec));
    fn.fold_an_mouse  = fullfile(local_disk,'an_data',sprintf('mouse_%s',mouse));
    fn.fold_an_sess   = fullfile(fn.fold_an_mouse,sprintf('sess_%s',sess));
    fn.basename_an    = sprintf('%s_%s_%s_',mouse,sess,rec);
    fn.exp_spikes   = fullfile(fn.fold_exp_data,   sprintf('%s_%s_%s_spikes.mat', mouse, sess, rec));
    fn.exp_trial    = fullfile(fn.fold_exp_data,   sprintf('%s_%s_%s_trial.mat', mouse, sess, rec));
end

%if station was specified
if nargin >3
    if isnumeric(stat)
        stat = sprintf('%02d', stat);
    end
    
    fn.fold_sd_data  = fullfile(stat_disk, sprintf('stat_%s',stat),'raw_data');
    fn.fold_sd_mouse = fullfile(fn.fold_sd_data,   sprintf('mouse_%s',mouse));
    fn.fold_sd_sess  = fullfile(fn.fold_sd_mouse,  sprintf('sess_%s',sess));
end

end
