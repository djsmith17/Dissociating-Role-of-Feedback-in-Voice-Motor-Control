function dfRunSubjAnalysis()
%Analyses the microphone data from the somatosensory perturbation
%experiment. Measures the change in f0 over each trial, and each run for a
%given participant. At the end it approximates a general response to
%inflation to be used in the auditory perturbation experiment

%Require the Signal Processing Toolbox

clear all; close all; clc
AVar.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
AVar.participants  = {'PureTone200'}; %List of multiple participants.
AVar.numPart       = length(AVar.participants);
AVar.runs          = {'AF6'};
AVar.numRuns       = length(AVar.runs);
AVar.debug         = 0;

dirs               = dfDirs(AVar.project);

for i = 1:AVar.numPart
    for j = 1:AVar.numRuns
        participant = AVar.participants{i};
        run         = AVar.runs{j};
        
        dirs.baselineData  = fullfile(dirs.RecData, participant, 'GT1', [participant 'GT1' 'DRF.mat']); %Where to find data
        dirs.SavFileDir    = fullfile(dirs.RecData, participant, run, [participant run 'DRF.mat']); %Where to find data
        dirs.SavResultsDir = fullfile(dirs.Results, participant, run); %Where to save full analyzed results
        
        dirs.InflaVarDir  = fullfile(dirs.RecData, participant, 'IV1');

        if exist(dirs.baselineData, 'file') == 0
            disp('ERROR')
            return
        end
        
        if exist(dirs.InflaVarDir, 'dir') == 0
            mkdir(dirs.InflaVarDir)
        end
        
        if exist(dirs.SavResultsDir, 'dir') == 0
            mkdir(dirs.SavResultsDir)
        end

        fprintf('Loading Files for %s %s\n', participant, run)
        load(dirs.baselineData) % Expect GT
        load(dirs.SavFileDir)   % Expect DRF
        
        AVar.expType = DRF.expParam.expType;
        [pF, iRF] = checkDRFExpType(AVar.expType);
        AudFlag = 1;
        
        %Initialize these so I can stop worrying about it
        niAn = []; niRes = [];
        auAn = []; auRes = [];
                
        bTf0b = GT.subjf0;
        [niAn, niRes] = dfAnalysisNIDAQ(dirs, DRF.expParam, DRF.DAQin, bTf0b, AudFlag, pF, iRF);
%         [auAn, auRes] = dfAnalysisAudapter(dirs, DRF.expParam, DRF.rawData, bTf0b, 1);

        
        dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [participant run 'ResultsDRF.mat']);
        if AVar.debug == 0
            fprintf('Saving Results for %s %s\n', participant, run)
            save(dirs.SavResultsFile, 'auAn', 'auRes', 'niAn', 'niRes')
        end
        
        if iRF == 1
            saveInflationResponse(dirs, niRes, participant, run, AVar.debug)
        end
    end
end
end

function [pF, iRF] = checkDRFExpType(expType)
% This sets different variables so that the analyzes are done differently. 

if strcmp(expType, 'Somatosensory Perturbation_Perceptual') == 1
    pF  = 1;      %Pressure Analysis Flag
    iRF = 1;      %Inflation Response Flag
else
    pF  = 0;      %Pressure Analysis Flag
    iRF = 0;      %Inflation Response Flag
end
end

function saveInflationResponse(dirs, res, participant, run, debug)

InflaVar = res.InflaStimVar;

dirs.InflaVarFile = fullfile(dirs.InflaVarDir, [participant 'IV1' 'DRF.mat']);
if debug == 0
    fprintf('Saving Inflation Stimulus Variables for %s %s\n', participant, run)
    save(dirs.InflaVarFile, 'InflaVar');
end
end