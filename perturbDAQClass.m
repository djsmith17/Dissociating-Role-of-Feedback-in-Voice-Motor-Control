classdef perturbDAQClass < handle
    %saveNIDAQ is a class that saves data brought in from the NI card. 
    
    properties
        curTrial
        recTime
        recData
        PerturbSigs
        nVS
    end
    
    methods
        function obj = perturbDAQClass(PerturbSigs, nVS)
            obj.curTrial = 1;
            obj.recTime  = [];
            obj.recData  = [];
            obj.PerturbSigs = PerturbSigs;
            obj.nVS      = nVS;
        end 
        
        function queuePertOutputData(obj, s)
            thisTrial = obj.curTrial;
            PertSig   = obj.PerturbSigs(:, thisTrial);
            NIDAQsig  = [PertSig obj.nVS];
            fprintf('Data for Trial %d cued\n', thisTrial)
            s.queueOutputData(NIDAQsig)
        end
        
        function updateNIDAQdata(obj, nTime, nData)
            obj.recTime = cat(1, obj.recTime, nTime);
            obj.recData = cat(1, obj.recData, nData);
        end
    end    
end