classdef saveNIDAQ < handle
    %saveNIDAQ is a class that saves data brought in from the NI card. 
    
    properties
        Time
        Data
    end
    
    methods
        function obj = saveNIDAQ()
            obj.Time = [];
            obj.Data = [];
        end        
        
        function updateNIDAQdata(obj, nTime, nData)
            obj.Time = cat(1, obj.Time, nTime);
            obj.Data = cat(1, obj.Data, nData);
        end
    end    
end

