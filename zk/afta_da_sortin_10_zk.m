%
% Data preparation after spike sorting with KlustaKwik
% zk, Dima Rinberg Lab, September 2013
%Use:
%   afta_da_sorting_09_zk(mouse,sess,serv)
%   mouse : mouse_id (string)
%   sess  : sess_id  (string)
%   serv  : server   (string) - optional, data server to push the data
%Generates the unit structures for a session, after the .clu files have
%been generated using ndm and KlustaKwik scripts kkrunlist.sh.
%Outupt is created in pr_data folder of the computational server.
%
% 2013-09-09
% it has been tested, it works ok with single and multiple records.
% one sticky point: if there are channels of the raw data that are skipped in the
% final data before the last channel corresponding to a site, then the
% sites are not going to be right.
% to fix this, i have to re-read the channel config, instead of extracting
% the channel names from the info file in function get_units(igroup).


function afta_da_sortin_10_zk(mouse, sess, serv)
global as

as.read_clu_info = @read_read_clu_info;
as.make_sp_files = @make_spk_files;
as.read_xml_file = @read_xml_file;
as.push_data     = @push_data;

if ~exist('mouse', 'var')
    mouse = 'ZKawakeM72nosniff';
    sess  = '001';
end

if ~exist('serv','var');
    serv='dataserver';
end

if nargin>2 && ~isempty(serv)
    push_spikes(mouse,sess,serv);
end

read_clu_info(mouse,sess);
make_spk_files(mouse,sess);
%for debugging
%read_xml_file(mouse, sess);
%get_spikes(mouse,sess,rec);
%TODO
%push_data(mouse,sess);
end


function make_spk_files(mouse,sess)
%reads the info file, to extact data on records
    fn = file_names(mouse, sess);
    %load(fn.ss_sess_info); %loads info structure
    
    clInfoBase=read_clu_info(mouse,sess);
    clInfoNew=clInfoBase;
    clInfoNew=rmfield(clInfoNew,'rec');

    
    for irec=1:numel(clInfoBase.rec)
        recname=clInfoBase.rec(irec).name
        %recname=info.rec(irec).name
        clInfoNew.rec(irec)=get_spikes(mouse,sess,recname,clInfoBase.rec(irec));
    end
    
    fprintf('Writing cluster info file %s\n',fn.clInfo_file);
    save(fn.clInfo_file,'-struct','clInfoNew');
    fprintf('>>> Done  making spike mat files for mouse %s session %s \n\n',mouse,sess);
    
end

function clInfoRec=get_spikes(mouse,sess,rec,clInfoRec)
% reads:
% 1 - xml file
% 2 - clu files
% 3 - res files
% 4 - spk files
% fills spike structure for ever spike in every cluster (got from xml file)
% sorts the spikes by event
fprintf('** Get spikes for rec %s: \n',rec);

fn = file_names(mouse, sess, rec)

if nargin<4 || isempty(clInfoRec)
    
    if ~exist(fn.clInfo_file, 'file')
        clInfo=read_clu_info(mouse,sess);
    else
        clInfo=load(fn.clInfo_file)
    end
    
    %%% this i want in a function get_spike_info
    %sessInfo=load(fn.ss_sess_info);
    if isnumeric(rec)
        rec = sprintf('%02d', rec);
    end
    
    indexOfRec=find(strcmpi(rec,{clInfo.rec.name}));
    clInfoRec=clInfo.rec(indexOfRec);
end

unit=struct();

fprintf('\n *** Getting cluster files for rec %s ***\n',rec);

clInfoRec.unitCount=0; %global counter of units (the final count of units is group independent);

for ig=1:clInfoRec.nGroup
    if isfield(clInfoRec,'groupInfo')
        iGroup=clInfoRec.sortGroups(ig);
    else
        iGroup=ig;
    end
    load_spikes(iGroup)
    get_sunits(iGroup)
end
clInfoRec.unit=unit;
clInfoRec.unitFn=fn.spikes;
fprintf('Writing units file %s\n',fn.spikes);
save(fn.spikes,'unit')
fprintf('--- Done getting spikes for rec %s\n\n',rec);
%save(['dbgspkinfo' rec '.mat'], 'clInfoRec')

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function load_spikes(igroup)
        %load all the spike files for a group into the clInfo struct.
        clInfoRec.group(igroup).nChan=length(clInfoRec.group(igroup).chan);        
        fprintf('Getting files for group %d/%d\n',igroup,clInfoRec.nGroup);
        cluFn=fullfile(fn.fold_ss_rec,['rec_' rec '.clu.' num2str(igroup)]);
        qltFn=fullfile(fn.fold_ss_rec,['rec_' rec '.qlt.' num2str(igroup)]);
        resFn=fullfile(fn.fold_ss_rec,['rec_' rec '.res.' num2str(igroup)]);
        spkFn=fullfile(fn.fold_ss_rec,['rec_' rec '.spk.' num2str(igroup)]);
        
        if ~exist(cluFn, 'file')
            error('no clu file')
        end
        clu=load(cluFn);
        clInfoRec.group(igroup).nClu=clu(1);
        clInfoRec.group(igroup).clu=clu(2:end);
        
        if ~exist(qltFn, 'file')
            defaultMask=3;
            fprintf('No spike description file for group %i:\n',igroup);
            fprintf('\t %s \n',qltFn);
            fprintf('\t ** Assuming all non-deleted spikes are good single units\n');
            qlt=defaultMask*ones(clu(1),1);
        else
            qlt=load(qltFn)
            clu(1)
            if isempty(qlt) || length(qlt)~=clInfoRec.group(igroup).nClu;
                error(['Spike quality description file does not match number of clusters for group ' num2str(igroup)])
            end
       end
        disp(qlt)
        clInfoRec.group(igroup).qltMasks=qlt;
        
        if ~exist(resFn, 'file')
            error('no res file')
        end
        res=load(resFn);
        clInfoRec.group(igroup).res=res;
        
        if ~exist(spkFn, 'file')
            error('no spk file')
        end
        spk=get_spk(spkFn);
        clInfoRec.group(igroup).spks=spk;
        
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        function spk=get_spk(spkFn)
            %it is a binary file with the same number of channels as the group
            if ~exist(spkFn, 'file')
                error('no spk file')
            end
            fid=fopen(spkFn);
            spk=fread(fid,[clInfoRec.group(igroup).nChan,inf],'int16');
            spk=spk';
        end %get_spk(spkFn)
        
    end %load_spikes(igroup)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    
    function get_units(igroup)
        %gets the units
        orderedUnits=unique(clInfoRec.group(igroup).clu);
        qltProps={'lightral','singleUnit','good'}; %properties, in the order of signifiance in the masking
        %i.e: {'lightral','multiunit; ,'good'} is masked 101=5
        % dont store artifacts and noise clusters (as maked by 0,1 in clu file)
        clInfoRec.group(igroup).qltMasks(orderedUnits<2)=[];
        orderedUnits(orderedUnits<2)=[];
        for ou=1:length(orderedUnits)
            cluId=orderedUnits(ou);
            clInfoRec.unitCount=clInfoRec.unitCount+1;
            cnt=clInfoRec.unitCount;
            unit(cnt).group=igroup;
            unit(cnt).clu=cluId;
            indexOfStamps= clInfoRec.group(igroup).clu==cluId;
            unit(cnt).nspikes=sum(indexOfStamps);
            unit(cnt).pkStamps=clInfoRec.group(igroup).res(indexOfStamps);
            %unit(cnt).stamps=unit(cnt).pkStamps-clInfoRec.group(igroup).peakSampleIndex;
            unit(cnt).times=unit(cnt).pkStamps/clInfoRec.sampling_freq*1000.;
            unit(cnt).nChans=length(clInfoRec.group(igroup).chan);
            unit(cnt).chans=clInfoRec.group(igroup).chan;
            % get the sites for those channels
            if ~isfield(clInfoRec,'groupInfo')
                for ich=1:length(unit(cnt).chans)
                    sites{ich}=clInfoRec.chanMap([clInfoRec.chanMap.num_rd]==unit(cnt).chans(ich)+1).name;
                end
            else
                sites={clInfoRec.group(igroup).chansDescr.name};
            end
        
            unit(cnt).sites=sites;
                       
            %get the avergage waveform for the spike
            spikes=clInfoRec.group(igroup).spks;
            iSpikes=find(indexOfStamps); %indexes of the spikes for this unit in the spk file
            iSpikes((iSpikes*32)>length(spikes))=[];%trim edge
            unit(cnt).nPoints=clInfoRec.group(igroup).nSamples;
            unitcuts=transpose(bsxfun(@plus,(iSpikes-1)*unit(cnt).nPoints,(1:unit(cnt).nPoints)));
            spikeWave=zeros(unit(cnt).nPoints,unit(cnt).nChans);
            %Debugging
%             unit(cnt).indexOfStamps=indexOfStamps;
            %unit(cnt).unitcuts=unitcuts;
            unict(cnt).ounits=spikeWave;
            %
            for ichan=1:unit(cnt).nChans
                spkSel=clInfoRec.group(igroup).spks(unitcuts,ichan);
                %unit(cnt).spkSel=spkSel;
                spkAvg=reshape(spkSel,unit(cnt).nPoints,[]);
                spkAvg=mean(spkAvg,2);
                spikeWave(:,ichan)=spkAvg;
            end
            
            unit(cnt).waveForm=spikeWave;
            %unit(cnt).oUnits=orderedUnits;
            %unit(cnt).qltMasks=clInfoRec.group(igroup).qltMasks;
            unit(cnt).qltMask=clInfoRec.group(igroup).qltMasks(orderedUnits==cluId);
            if isnumeric(unit(cnt).qltMask)
                unmasked=de2bi(unit(cnt).qltMask,'left-msb',numel(qltProps));
                for im=1:numel(qltProps)
                    unit(cnt).(qltProps{im})=unmasked(im);
                end
            end
            unit(cnt).sorted=(cluId>1);
        end
    end %get_units(spkFn)
end

function clInfo=read_clu_info(mouse,sess)
% reads session info and gets the clu info for all the recs in the session
fn = file_names(mouse, sess);

fprintf('Getting info file %s\n',fn.ss_sess_info);
load(fn.sess_info);
clInfo=struct();

fprintf('* Will get spikes for mouse %s, session %s ...\n',mouse,sess);
fprintf('(*----------------------------------------------------*) \n');

for irec=1:numel(info.rec)
    fprintf('Rec %s : \n',info.rec(irec).name);
    fn = file_names(mouse, sess,info.rec(irec).name);
    clInfo.mouse=info.mouse; %for later consistency checks
    clInfo.sess=info.sess; %for later consistency checks
    %if there was channel grouping information in the info file, use that
    %otherwise, read it from the xml file
    clInfo.rec(irec).xmlRead=true;
    xmlInfo=read_xml_file(fn.ss_xml);
    xmFiN=fieldnames(xmlInfo);
    for ifield=1:numel(xmFiN)
        clInfo.rec(irec).(xmFiN{ifield})=getfield(xmlInfo,xmFiN{ifield});
    end
    
    %the newer versions of post_exp_processing have a structure of channel
    %grouping information.
    %if such structure exists, take it from there.
    if isfield(info.rec,'chGroupInfo')
        fprintf('Found spike group information in the session info file\n')
        clInfo.rec(irec).groupInfo=true;
        clInfo.rec(irec).nGroup=numel(find([info.rec(irec).chGroupInfo.sortable]));
        clInfo.rec(irec).sortGroups=find([info.rec(irec).chGroupInfo.sortable]);
        for isg=1:clInfo.rec(irec).nGroup
            groupIdx=clInfo.rec(irec).sortGroups(isg);
            clInfo.rec(irec).group(isg).groupNumber = groupIdx;
            clInfo.rec(irec).group(isg).name        = info.rec(irec).chGroupInfo(groupIdx).name;
            clInfo.rec(irec).group(isg).chan        = info.rec(irec).chGroupInfo(groupIdx).chan;
            clInfo.rec(irec).group(isg).nChan       = info.rec(irec).chGroupInfo(groupIdx).nCh;
            clInfo.rec(irec).group(isg).chansDescr  = info.rec(irec).chan(clInfo.rec(irec).sortGroups);
        end
    end
    
    clInfo.rec(irec).name=info.rec(irec).name;
    clInfo.rec(irec).sampling_freq=info.rec(irec).sampling_freq;
    clInfo.rec(irec).chanMap=info.rec(irec).chan;
end
save(fn.clInfo_file,'-struct','clInfo');
fprintf('\n--- Clu info for mouse %s, session %s read.\n',mouse,sess);
end

function [clInfo] = read_xml_file(xml_file_path)
%get spike groups from the xml file
%read sequentially, if <spikeDetection> key
fprintf('  *** Reading xml file %s ***\n',xml_file_path);
clInfo   = struct();
irec   = 0;
if ~exist(xml_file_path, 'file')
    error('no xml file')
end
fid = fopen(xml_file_path, 'r');

neof   = true;
plevel = 0;
addToStru=0; %flag for adding the
while neof
    % read a line from xml file
    tl = fgetl(fid);
    % if end of the file, finish the prorgam,
    % if the line starts with ':' - read the parameter
    if tl== -1 | strcmp(strtrim(tl),'OPTIONAL PARAMETERS') | strcmp(strtrim(tl),'</spikeDetection>')
        neof =0;
    elseif ~isempty(tl) && ~isempty(strtrim(tl))
        tl=strtrim(tl);
        if tl(1)=='<' %it is a valid xml line
            tl=tl(2:end);
            [key,rol]=strtok(tl,'>');
            rol(end)='';
            %fprintf('key %s,  rol %s \n',key,rol);
            if ~isempty(rol)
                % key is a key to to a value get the value.
                value=strtok(rol(2:end),'<');
                %fprintf('**** %s  %s, plevel %02d\n',key,value,plevel);
                if ~isnan(value)
                    value = str2double(value);
                end
            elseif ~(key(end)=='/')
                %key is a key to a level jump
                if key(1)=='/'
                    plevel=plevel-1;
                else
                    plevel=plevel+1;
                end
                %fprintf('---- key %s,  plevel %02d \n',key,plevel);
            end
            
            %depending on what the key is, add it to the structure
            %spikeDetection is level 3
            %channelGroups is level 4
            %group is level 5
            %channels is level 6
            %nSamples is level 6
            %peakSampleIndex is level 6
            %nFeatures is level 6
            %is spikeDetection is turned on, add the following
            %fields to the structure
            if strcmp(key,'spikeDetection')
                addToStru=1;
                pl0=plevel;
                igroup=0;
            elseif strcmp(key,'/spikeDetection')
                addToStru=0;
                disp(clInfo);
            end
            
            if addToStru==1 && ~(key(1)=='/');
                switch key
                    case 'group'
                        igroup=igroup+1;
                        value=igroup;
                        key='groupNumber';
                        fprintf('\t Group %d: \n',igroup);
                    case 'channel'
                        ich=ich+1;
                end
                
                plrel=plevel-pl0;
                %fprintf('---- key %s,  plevel %02d, plrel %02d \n',key,plevel,plrel);
                switch plrel
                    case 2
                        clInfo.group(igroup).(key)=value;
                        ich=0;
                    case 3
                        if(strcmp(key,'channel'))
                            clInfo.group(igroup).chan(ich)=value;
                            fprintf('\t \t Channel %02d : %02d \n',ich ,value);
                        end %switch
                end %if addtostru=1
            end % if addToStru==1 && ~(key(1)=='/');
        end %if tl(1)=='>'
    end %if tl(1)== -1
end  %while
fclose(fid);
clInfo.nGroup=numel(clInfo.group);
end
