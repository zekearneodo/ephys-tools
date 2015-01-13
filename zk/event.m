%a parent class for events
classdef event < handle
    properties
        mouse = '';
        sess  = '';
        rec   = '';
        %
        rawFile  = '';
        chanName = '';
        name     = '';
        samplingFrequency = 19531;
        %
        eventTimes = [];
    end

    methods (Abstract = true)
        function obj = event(obj, mouse, sess, rec, name, varargin)
            inPar=inputParser;
            inPar.addRequired('mouse')
            inPar.addRequired('sess')
            inPar.addRequired('rec',@ischar)
            inPar.addRequired('name',@ischar);
            
            inPar.addParamValue('chanName','',@(x) ischar(x));
            inPar.addParamValue('name','',@(x) ischar(x));
            inPar.addParamValue('samplingFrequency',19531,@(x) isscalar(x))
            
            inPar.parse(mouse,sess,rec,name,varargin{:});
            obj.mouse = inPar.Results.mouse;
            obj.sess  = inPar.Results.sess;
            obj.rec   = inPar.Resuts.rec;
            obj.name  = inPar.Resutls.name:
            
            if isempty(inPar.Results.chanName)
                obj.chanName = inPar.Results.Name;
            else
                obj.chanName = obj.name;
            end
            
            obj.samplingFrequency = inPar.Results.samplingFrequency;
            
        end
       
       %function load_stream = @load_raw_data(obj)
       %function 
    end
end