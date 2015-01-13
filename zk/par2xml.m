function parentNode=par2xml(xo,stru,parentNode)
% gets a parameter structure and turns it into an xml node with the right
% format for the ndm xml files, parameter programs section.
% xo   = root DOM of an xml file.
% stru = (structure) ndm-ready parameter structure.
% parentNode = a parent DOM node to which the created node will be appended
% as a childNode. if no parentNode is entered, it creates an empty
% parentNode and recursively appends the subnodes as it goes through the
% structure.
% example:
% extractSpikesNode=par2xml(xmlread(fn.ss_xml),par_ndm_extractspikes)
% (creates a node with the parameters for program ndm_extractspikes)
%par2xml
%go through the parameter structure,
%create the node with Name=key, value=value.
%if the struct has children, do the same for every child item
keys=fields(stru);
if nargin<3 || isempty(parentNode)
    parentNode=xo.createElement(keys{1});
    stru=stru.(keys{1});
    if numel(keys)>1
        warning('Root of parameter structure has multiple entries. It should have only one ("parameter")')
    end
    parentNode=par2xml(xo,stru,parentNode);
    %parentNode.appendChild(newNode);
    return 
end

for ik=1:numel(keys)
    %make node with name of par
    subStru=stru.(keys{ik});
    newNode=xo.createElement(keys{ik});
    
    %if it has a value, set the value (node name=key, node
    %TextContent=value)
    if ~isstruct(subStru)
        value=stru.(keys{ik});
        if isnumeric(value)
            value=num2str(value);
        end
        newNode.setTextContent(value);
    else
        %if it has a structure, make a node out of that structure
        for iss=1:numel(subStru)
            newNode=par2xml(xo,subStru(iss),newNode);
        end
    end
    parentNode.appendChild(newNode);
end

end