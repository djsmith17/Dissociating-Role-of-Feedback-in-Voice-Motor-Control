function dfRunSubjAnalysis()
%Analyses the microphone data from the somatosensory perturbation
%experiment. Measures the change in f0 over each trial, and each run for a
%given participant. At the end it approximates a general response to
%inflation to be used in the auditory perturbation experiment

%Require the Signal Processing Toolbox

clear all; close all; clc
AVar.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
AVar.participants  = {'Pilot24', 'Pilot25', 'Pilot26', 'Pilot22'}; %List of multiple participants.
AVar.numPart      = length(AVar.participants);
AVar.runs         = {'SF1'; 'SF2'; 'SF3'; 'SF4'};
AVar.numRuns      = length(AVar.runs);
AVar.debug        = 0;

dirs               = dfDirs(AVar.project);

for i = 1:AVar.numPart
    for j = 1:AVar.numRuns
        participant = AVar.participants{i};
        run         = AVar.runs{j};
        
        dirs.baselineData  = fullfile(dirs.RecData, participant, 'GT1', [participant 'GT1' 'DRF.mat']); %Where to find data
        dirs.SavFileDir    = fullfile(dirs.RecData, participant, run, [participant run 'DRF.mat']); %Where to find data
        dirs.SavResultsDir = fullfile(dirs.Results, participant, run); %Where to save results

        if exist(dirs.baselineData, 'file') == 0
            disp('ERROR')
            return
        end
        
        if exist(dirs.SavResultsDir, 'dir') == 0
            mkdir(dirs.SavResultsDir)
        end

        fprintf('Loading Files for %s %s\n', participant, run)
        load(dirs.baselineData)
        load(dirs.SavFileDir)
        
        bTf0b = GT.subjf0;
        [niAn, niRes] = dfAnalysisNIDAQ(dirs, DRF.expParam, DRF.DAQin, bTf0b, 1);
        [auAn, auRes] = dfAnalysisAudapter(DRF.expParam, DRF.rawData, DRF.DAQin);

        dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [participant run 'ResultsDRF.mat']);
        if AVar.debug == 0
            fprintf('Saving Results for %s %s\n\n', participant, run)
            save(dirs.SavResultsFile, 'auAn', 'auRes', 'niAn', 'niRes')
        end
    end
end
end