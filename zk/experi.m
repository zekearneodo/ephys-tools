classdef parameter < handle
    properties
        name  = '';
        value = '';
    end
    
    methods
        function obj = parameter(name, value)
            obj.name  = name;
            obj.value = value;
        end
        
        function obj = print(fid)
            if nargin>0 && ~isempty(fid)
                fprintf(fid,'%s = %d\n',obj.name,obj.value);
            else
                fprintf('%s = %d\n',obj.name,obj.value)
            end
            
        end
    end
    
end