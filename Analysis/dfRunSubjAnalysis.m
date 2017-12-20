function dfRunSubjAnalysis()
%Analyses the microphone data from the somatosensory perturbation
%experiment. Measures the change in f0 over each trial, and each run for a
%given participant. At the end it approximates a general response to
%inflation to be used in the auditory perturbation experiment

%Require the Signal Processing Toolbox

clear all; close all; clc
AVar.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
AVar.participants  = {'Pilot24'}; %List of multiple participants.
AVar.numPart       = length(AVar.participants);
AVar.runs          = {'SF1'};
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
        load(dirs.baselineData)
        load(dirs.SavFileDir)
        
        bTf0b = GT.subjf0;
        [niAn, niRes] = dfAnalysisNIDAQ(dirs, DRF.expParam, DRF.DAQin, bTf0b, 1);
        [auAn, auRes] = dfAnalysisAudapter(dirs, DRF.expParam, DRF.rawData, bTf0b, 1);

        InflaVar = niRes.InflaStimVar;

        dirs.InflaVarFile   = fullfile(dirs.InflaVarDir, [participant 'IV1' 'DRF.mat']);
        dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [participant run 'ResultsDRF.mat']);
        if AVar.debug == 0
            fprintf('Saving Results for %s %s\n', participant, run)
            save(dirs.SavResultsFile, 'auAn', 'auRes', 'niAn', 'niRes')
            fprintf('Saving Inflation Stimulus Variables for %s %s\n', participant, run)
            save(dirs.InflaVarFile, 'InflaVar');
        end
    end
end
end