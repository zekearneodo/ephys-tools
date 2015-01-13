function sniff_correct_007(sess,rec)

% check that the file is not one that has been corrected already! 
switch sess
    case '005'
        done = {'a'};
    case '006'
        done = {'a'}
    case '007'
        done = {'a','b','c','d'}
end
if strcmp(rec,done)
    error('Hold on!! KPawakeM72 sess %s rec %s was already corrected!',sess,rec);
end


includePath=fullfile(fileparts(pwd),'current','include');
addpath(includePath);

dm = data_management_tools_032();
fn = dm.file_names('KPawakeM72', sess, rec);

load(fn.trial)

for ikt = 1:length(trial)
    zeros = trial(ikt).sniffZeroTimes;
    parabs = trial(ikt).sniffParabZeroTimes;
    TTL = trial(ikt).sniffTTLTimes;
    
    trial(ikt).sniffZeroTimes = -12 + zeros;
    trial(ikt).sniffParabZeroTimes = -12 + parabs;
    trial(ikt).sniffTTLTimes = -12 + TTL;
    
end
    save(fn.trial,'trial')

end