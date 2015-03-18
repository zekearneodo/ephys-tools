%Testing of the selialNumbers class and the virtualserial classes
% ZK
% December 2014

clear serialNumbers serialdata virtualserial
close all
%get me a serial from a mouse
mouse = 'ZKawakeM72';
sess  = '007';
rec   = 'c';
irun   = 1;
chanName = 'trNum';

%make the serial data chunk
serialNumbers = serialdata('serial',{300,8,0,1,1000})

fn=file_names(mouse,sess,rec)
q=load(fn.ss_sess_info)
sInfo=q.info
recInfo=sInfo.rec(find(strcmpi(rec,{sInfo.rec.name})))
runInfo = recInfo.run([recInfo.run.num]==irun)
runRawFile=fullfile(fn.fold_rd_sess,runInfo.ephys_data)

trnChan=recInfo.chan(strcmpi(chanName,{recInfo.chan.name}))

fprintf('Reading serial trial numbers channel (%s)...',chanName);
serialNumbers.read_stream_bin(runRawFile,trnChan.num_rd,recInfo.nChan_rd,'int16')
fprintf('done.\n');

serialNumbers.clip_stream()
%plot the things
figure
hold on
plot(serialNumbers.stream(1:200000))
ylim([-1.5 1.5])