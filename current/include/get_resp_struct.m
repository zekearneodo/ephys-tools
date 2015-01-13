function resp = get_resp_struct(unit_filename,stimType)
% fn=file_names();

% % add option to get filepath if input arg is cellId
% cellFilePath=fullfile(fn.fold_unit_db, sprintf('%s.mat',cellId));

unitInfo=load(unit_filename);

fn=file_names(unitInfo.mouse,unitInfo.sess,unitInfo.rec);
cluList='';
for iu=1:numel(unitInfo.clu)
    cluList=[cluList num2str(unitInfo.clu(iu),'%02d')];
end

unitRespFn     = sprintf('%s%s_units%s_resp.mat',fn.basename_an,stimType,cluList);
unitRespFullFn = fullfile(fn.fold_an_sess,unitRespFn);
% unitResp       = struct;

load(unitRespFullFn);


end