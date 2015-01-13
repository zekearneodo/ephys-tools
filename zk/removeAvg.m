%very basic script to remove the first candidates of common noise
%step 0:
%find the spike candidates that appear exactly at the same time in all the
%groups, and make them a new cluster.
%step 1:
%instead of removing them from the pool, substract the average in all the
%ephys channels; then run the detection again
%step 2:
%get the raw data file, find the PCs for the groups at the spike
%candidates. If the first PC is the same for all the groups, substract the
%projection to each individual channel.

%Begin with step 1:
%if this works, it will be done in the ss_prep step
%includes:
ft=file_tools();

%Open the fet file:
mouse='SDbehavingM72';
sess=1;
rec='a';
irec=1;

fn=file_names(mouse,sess,rec);
q=load(fn.ss_sess_info);
thisRec=q.info.rec(irec);

%list all the channels that are in any group
totalChans=thisRec.nChan;
ephysGroups=thisRec.chGroupInfo(find([thisRec.chGroupInfo.sortable]));
ephysChans=[ephysGroups.chan]+1;

%open the file:
recFn=fullfile(fn.fold_ss_rec,sprintf('rec_%s.fil',rec))
dFid=fopen(recFn,'r');
%read one by one the ephys channels and compute the avg channel;
fprintf('++ Reading all ephys channels and computing average...\n');
fprintf('   ');
avgCh=ft.read_analog_channel(dFid,ephysChans(1),totalChans);
avgAux=zeros(size(avgCh));
for iCh=ephysChans(2:end)
    fprintf('   ');
    avgAux=ft.read_analog_channel(dFid,iCh,totalChans);
    avgCh=avgCh+avgAux;
end
fclose(dFid);
avgCh=avgCh/numel(ephysChans);
fprintf('-- %d channels read; avg computed.\n',numel(ephysChans));

%now save every channel, with the avg subtracted.
fprintf('xx Reading all ephys channels, subtracting average and saving in rec_%s.cle...\n',rec);
wFid=fopen(fullfile(fn.fold_ss_rec,sprintf('rec_%s.cle',rec)),'w');
dFid=fopen(recFn,'r');
for iCh=ephysChans
    fprintf('   ');
    avgAux=ft.read_analog_channel(dFid,iCh,totalChans);
    avgAux=avgAux-avgCh;
    fprintf('   ');
    ft.write_analog_channel(avgAux,wFid,0,iCh,totalChans);
end

fclose(dFid);
fclose(wFid);
fprintf('-- %d channels read; avg subrtracted and saved.\n',numel(ephysChans));