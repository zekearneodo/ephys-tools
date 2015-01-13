function [baseFolder]=get_local_disk()
%sets up base folder for the experiment
%if it is lunux or mac, checks the hostname
persistent local_disk;

if isempty(local_disk)
%sets up base folder for the experiment
%if it is lunux or mac, checks the hostname
[~,computerName]=system('hostname');
%if its the computational server, load the configurations.
%careful when you replace this, keep in mind it will affect all the
%computers.
%Here, of all places, changes MUST be backwards compatible
fprintf('Setting base folder for computer is %s)\n',computerName);
switch strtrim(computerName)
    case 'flipper'
        %the file server
        config_fold = '/usr/local/kluster/config';
        cf=fopen(fullfile(config_fold,'experiment_folder'),'r');
        local_disk=fgetl(cf);
        fclose(cf);
        %local_disk='/home/zeke/data';
        fn.xml_template=fullfile(config_fold,'template.xml');
        stat_disk='/stations';
    case 'Dima-MacBookPro4.local'
        %zk notebook
        local_disk  = '/Users/zeke/experiment';
    case 'kristina-macpro.ersp.nyumc.org';
        local_disk = '/Volumes/spikefolder';
    case 'kristinikissmbp.vz30.nyumc.org';
        local_disk = '/Users/kpenikis/Documents/Rinberg/Experiment';
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
baseFolder=local_disk;
end %get_local_disk()