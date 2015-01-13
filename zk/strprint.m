function meta = strprint( stru,meta)
% prints a ndm_structure in ndm xml format to a file.
% if the structure is parameters, it adds the ndm header <?xml
% version='1.0'?> in level 0,
% and pastes the standard, post-ephys parameters part of the xml file in
% the end.
%
% checks how many fields the structure has.
% for each of them: 
% if it has a val field, it prints it.
% if none, it prints the null key
if nargin<2
    level=0;
    struName=inputname(1);
    meta.fWid=fopen('ndm.xml','w');
    fprintf(meta.fWid,'<?xml version=\''1.0\''?>\n');
else
    level=meta.level+1;
    struName=meta.name;
end

%make the indent
levelTab='';
for it=1:level
    levelTab=[levelTab ' '];
end

%get the contents of the struct and print them
keys=fields(stru);
%if its a null value field, print it as a null value field
if numel(keys)==0
    fprintf(meta.fWid,'%s<%s/>\n',levelTab,struName);

else
    %if it is a non-null field, open the key
    %check if it has an option
    optStr='';
    if sum(strcmp('opt',keys))
        if ~isempty(stru.opt) 
            optStr=[' ' stru.opt];
        end
        keys(strcmp(keys,'opt'))=[];
    end
    
    fprintf(meta.fWid,'%s<%s%s>',levelTab,struName,optStr);    
    %if it is a value, print the value and get ready to close
    if sum(strcmp('val',keys))
        fprintf(meta.fWid,'%s',num2str(stru.val));
        levelTab='';
    else
        fprintf(meta.fWid,'\n');
        keys(strcmp(keys,'val'))=[];
        %if it has other fields, print them and print carriage return
        for ik=1:numel(keys)
            %fprintf('printkey %s<%s>\n',levelTab,keys{ik});
            key=keys{ik};
            subStruArray=stru.(key);
%             subStruArray(1)
            for ia=1:length(subStruArray)
                meta.name=keys{ik};
                meta.level=level;
                subStru=subStruArray(ia);
                meta=strprint(subStru,meta);
            end
                
        end
    end
    if ~strcmp(struName,'parameters')
        fprintf(meta.fWid,'%s</%s>\n',levelTab,struName);
    else
        %insert template of rest of options' file
        fRid=fopen('template.xml','r');
        fseek(fRid,0,-1);
        while ~feof(fRid)
            tline=fgetl(fRid);
            fprintf(meta.fWid,'%s\n',tline);
        end
        fclose(fRid);
        fclose(meta.fWid);
    end
end
%when dne goin through, decrease level in 1
level=level-1;

end

