classdef perturbDAQClass < handle
    %saveNIDAQ is a class that saves data brought in from the NI card. 
    
    properties
        curTrial
        recTime
        recData
        trialTime
        trialData
        PerturbSigs
        nVS
    end
    
    methods
        function obj = perturbDAQClass(PerturbSigs, nVS)
            obj.curTrial    = 1;
            obj.recTime     = [];
            obj.recData     = [];
            obj.trialTime   = [];
            obj.trialData   = [];
            obj.PerturbSigs = PerturbSigs;
            obj.nVS         = nVS;
        end 
        
        function queuePertOutputData(obj, s)
            if s.IsLogging == 0 && s.TriggersRemaining > 0
                thisTrial = obj.curTrial;
                PertSig   = obj.PerturbSigs(:, thisTrial);
                NIDAQsig  = [PertSig obj.nVS];
                fprintf('Data for Trial %d cued\n', thisTrial)
                s.queueOutputData(NIDAQsig)
            end
        end
        
        function updateNIDAQdata(obj, s, nTime, nData)
            fprintf('Received Data\n')
%             fprintf('IsLogging = %d, IsRunning = %d, IsDone = %d\n', s.IsLogging, s.IsRunning, s.IsDone)
            obj.recTime = cat(1, obj.recTime, nTime);
            obj.recData = cat(1, obj.recData, nData);
        end
        
        function storeTrialData(obj, s)
            if s.IsLogging == 0 && isempty(obj.recData) == 0
                fprintf('Logging Previous Trial\n')
                obj.trialTime = cat(2, obj.trialTime, obj.recTime);
                obj.trialData = cat(3, obj.trialData, obj.recData);
                obj.recTime = [];
                obj.recData = [];
            end
        end
    end    
end