function make_x_file(mouse, sess)
%make the structure for an experiment's hd5 file
%group structure of the file is
%sess
%  rec
%    run
%      events    
fn=file_names(mouse,sess);
q=load(fn.sess_info);
sInfo=q.info;

plist = 'H5P_DEFAULT';
%make the file
fcpl = H5P.create('H5P_FILE_CREATE');
fapl = H5P.create('H5P_FILE_ACCESS');
fid = H5F.create(fn.sess_h5,'H5F_ACC_TRUNC',fcpl,fapl);


for iRec=1:numel(sInfo.rec)
    %create the rec group name
    recGName=sInfo.rec.name;
    recGID = H5G.create(fid,recGName,plist,plist,plist);
    for iRun=1:numel(sInfo.rec(iRec).run)
        runGName=num2str(sInfo.rec(iRec).run(iRun).num,'%02d');
        runGID = H5G.create(recGID,runGName,plist,plist,plist);
        make_groups(runGID,{'streams','events'});
        H5G.close(runGID);
    end
    H5G.close(recGID);
end
H5F.close(fid);

    % ~~~~~~~~~~~~~~~~~~~~~~
    function make_groups(parentGid,groupList)
        for ig=1:numel(groupList)
            plist = 'H5P_DEFAULT';
            gid = H5G.create(parentGid,groupList{ig},plist,plist,plist);
            H5G.close(gid);
        end
    end

end
