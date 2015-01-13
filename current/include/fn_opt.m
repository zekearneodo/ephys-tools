%this function makes the file structure for one element of ephys data.
%receives parameters of the recording, and returns a structure with file
%and folder names.
%Oct 2013, DR
%Modified by ZK

%It should be placed in /basefolder/ephysDataManagement/current/include
%functions calling it should include that folder in the path by calling
%addpath('/basefolder/ephysDataManagement/current/include') 
%or something like that.

function fn = fn_opt(mouse, sess, rec, stat)
%this function makes the file structure for one element of ephys data.
%receives parameters of the recording, and returns a structure with file
%and folder names.

fn.local_disk=get_local_disk();

fn.fold_rd_data  = fullfile(fn.local_disk, 'raw_data');
fn.fold_ss_data  = fullfile(fn.local_disk, 'ss_data');
fn.fold_pr_data  = fullfile(fn.local_disk, 'pr_data');
fn.fold_config   = fullfile(fn.local_disk, 'SpikeGL_config');
if exist('stat_disk','var')
    fn.fold_sd_data  = fullfile(stat_disk, 'var');
end


if nargin > 0
    if isnumeric(mouse)
        mouse = sprintf('%04d', mouse);
    end
    if isnumeric(sess)
        sess = sprintf('%03d', sess);
    end
    fn.mouse   = mouse;
    fn.sess    = sess;
    fn.rec     = '';
    fn.fold_rd_mouse = fullfile(fn.fold_rd_data,  sprintf('mouse_%s', mouse));
    fn.fold_rd_sess  = fullfile(fn.fold_rd_mouse, sprintf('sess_%s', sess));
    fn.log           = fullfile(fn.fold_rd_sess,  sprintf('log_%s_%s.txt', mouse, sess));
    
    fn.fold_ss_sess  = fullfile(fn.fold_ss_data,  sprintf('ss_%s_%s', mouse, sess));
    fn.fold_pr_sess  = fullfile(fn.fold_pr_data,  sprintf('m%s_%s', mouse, sess));
    fn.ss_sess_info  = fullfile(fn.fold_ss_sess,  sprintf('%s_%s_sess_info.mat', mouse, sess));
    fn.sess_info     = fullfile(fn.fold_pr_sess,  sprintf('%s_%s_sess_info.mat', mouse, sess));
end

%if rec was specified
if nargin >2 && ~isempty(rec)
    if isnumeric(rec)
        rec = sprintf('%02d', rec);
    end
    
    fn.rec = rec;
    fn.fold_ss_rec  = fullfile(fn.fold_ss_sess,   sprintf('rec_%s',rec));
    fn.rd_rec_bn    = sprintf('%s',rec);
  
    fn.ss_rec       = fullfile(fn.fold_ss_rec,    sprintf('rec_%s.dat',rec));
    fn.ss_xml       = fullfile(fn.fold_ss_sess,   sprintf('rec_%s.xml',rec));
    fn.rsm_data     = fullfile(fn.fold_pr_sess,   sprintf('%s_%s_%s_rsm.mat', mouse, sess, rec));
    fn.trial        = fullfile(fn.fold_pr_sess,   sprintf('%s_%s_%s_trial.mat', mouse, sess, rec));
    fn.spikes       = fullfile(fn.fold_pr_sess,  sprintf('%s_%s_%s_spikes.mat', mouse, sess, rec));
    fn.sniffs       = fullfile(fn.fold_pr_sess,  sprintf('%s_%s_%s_sniff.mat', mouse, sess, rec));
    
    fn.fold_an_mouse  = fullfile(local_disk,'an_data',sprintf('mouse_%s',mouse));
    fn.fold_an_sess   = fullfile(fn.fold_an_mouse,sprintf('sess_%s',sess));
    fn.basename_an    = sprintf('%s_%s_%s_',mouse,sess,rec);
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
