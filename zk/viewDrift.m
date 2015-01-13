
mouse='ZKflowtest';
sF=19385;
subsample=19;
actaulFreq=round(sF/subsample);

chan=1;
totalChan=1;

sess=1;
rec='d';
run=1;

flowStart=900;
flowEnd=750;


%get the file, read it subsampled
fn=file_names(mouse,sess,rec);
binFileName=fullfile(fn.fold_rd_sess,sprintf('%s_%02d_data.bin',rec,run));
fid=fopen(binFileName,'r');
fseek(fid, 2*(chan-1), 'bof');
fprintf('Reading the flow channel in %s ...\n',binFileName);
%read subsampled (for instance to ~ms, approximately)
Y = fread(fid, inf, 'int16', (totalChan-1)*2+subsample*totalChan*2);
%With this subsample it is down to one sapmple per second;

%get minute-wise vector, and its range;
minutes=floor(length(Y)/60/1000)
Yminutes=reshape(Y(1:minutes*60000),60000,minutes);
fluctMinutes=[min(Yminutes,[],1)' max(Yminutes,[],1)'];
meanMinutes=mean(Yminutes,1);
sdMinutes=[(meanMinutes + sqrt(var(Yminutes,1)))' (meanMinutes - sqrt(var(Yminutes,1)))'];

%scale, visualize
b=(mean(meanMinutes(end-5:end))-mean(meanMinutes(1:5)))/(flowEnd-flowStart);
a=(mean(meanMinutes(end-5:end))/flowEnd-mean(meanMinutes(1:5))/flowStart)/(1/flowEnd-1/flowStart)
%a=mean(meanMinutes(end-5:end))-b*flowEnd;

%now scale to flow
meanMinutesFlow=(meanMinutes-a)/b;
fluctMinutesFlow=(fluctMinutes-a)/b;
sdMinutesFlow=(sdMinutes-a)/b;

%plot the result
flowFig=figure()
plot(meanMinutesFlow)
hold
plot([fluctMinutesFlow sdMinutesFlow])
title(sprintf('Flow fluctuations, run %03d-%s-%02d',sess,rec,run), 'FontSize', 10, 'FontWeight', 'bold')
if ~exist(fn.fold_pr_sess)
    mkdir(fn.fold_pr_sess)
end
print('-dpdf',fullfile(fn.fold_pr_sess,sprintf('%s_%02d_flow.pdf',rec,run)));
hold off

flowSpan=max(diff(fluctMinutesFlow'))
flowAvgSpan=mean(diff(fluctMinutesFlow'))