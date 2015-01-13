function [fileNames, recInfo] = get_rec_info(mouse,sess,rec)
%gets the info structure of a mouse, sess, rec
%input is
%mouse : string
%sess  : string/int
%rec   : string
recInfo = struct();
fileNames = file_names(mouse,sess,rec);
load(fileNames.ss_sess_info);

recIndex = find(strcmpi('a',{info.rec.name}));

if isempty(recIndex)
    warning('Did not find info for %s ' , fileNames.basename_an(1:end-1));
elseif numel(recIndex)>1
    warning('Found more than one match for info of %s, returning first match ' , fileNames.basename_an(1:end-1));
    recIndex = recIndex(1);
end
    
recInfo = info.rec(recIndex);

end