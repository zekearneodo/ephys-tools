classdef virtualserial < handle
%a virtual serial, to extract information sent by arduino in serial form
    properties(SetAccess = protected)
        baudRate = 19200;
        nBits    = 8;
        parity   = 0;
        stopBits = 1;
        timeOut  = 1000; % (timeout  between bytes in ms)
    end
    
    methods
        function obj = virtualserial(baudRate,nBits,parity,stopBits,timeOut)
            % the constructor
            if nargin>0
                obj.baudRate = baudRate;
                obj.nBits    = nBits;
                obj.parity   = parity;
                obj.stopBits = stopBits;
                obj.timeOut  = timeOut;
            end
        end
        
    end
        
end

