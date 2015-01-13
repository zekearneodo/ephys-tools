function out=execute_fun(arg1,arg2,varargin)

p=inputParser;

addRequired(p,'arg1',@ischar)
addRequired(p,'arg2',@ischar)


addParamValue(p,'f','file_names',@ischar)
addParamValue(p,'path','',@ischar)

disp(nargin)

if rem(nargin,2)~=0
    addOptional(p,'arg3','',@ischar)
    parse(p,arg1,arg2,varargin{:})
else
    parse(p,arg1,arg2,varargin{:})
end

disp(p)
out=p
% fh=str2func(f)
% out=fh(arg1,arg2);
end