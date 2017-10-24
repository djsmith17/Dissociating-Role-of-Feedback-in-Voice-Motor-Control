classdef perturbDAQClass < handle
    %saveNIDAQ is a class that saves data brought in from the NI card. 
    
    properties
        recTime
        recData
    end
    
    methods
        function obj = perturbDAQClass()
            obj.recTime = [];
            obj.recData = [];
        end        
        
        function updateNIDAQdata(obj, nTime, nData)
            obj.recTime = cat(1, obj.recTime, nTime);
            obj.recData = cat(1, obj.recData, nData);
        end
    end    
end