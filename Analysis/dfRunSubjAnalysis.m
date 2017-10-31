function dfRunSubjAnalysis()
%Analyses the microphone data from the somatosensory perturbation
%experiment. Measures the change in f0 over each trial, and each run for a
%given participant. At the end it approximates a general response to
%inflation to be used in the auditory perturbation experiment

%Require the Signal Processing Toolbox

clear all; close all; clc
AVar.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
AVar.participant  = 'Pilot0'; %List of multiple participants.
AVar.run          = 'SF3';
AVar.debug        = 1;

dirs               = dfDirs(AVar.project);
dirs.SavFileDir    = fullfile(dirs.SavData, AVar.participant, AVar.run, [AVar.participant AVar.run 'DRF.mat']); %Where to find data
dirs.SavResultsDir = fullfile(dirs.Results, AVar.participant, AVar.run); %Where to save results
dirs.InflaRespFile = fullfile(dirs.SavData, AVar.participant, [AVar.participant '_AveInflaResp.mat']);

if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

fprintf('Loading Files for %s %s\n', AVar.participant, AVar.run)
load(dirs.SavFileDir)

% [auAn, auRes] = dfAnalysisAudapter(DRF.expParam, DRF.rawData, DRF.DAQin);
[niAn, niRes] = dfAnalysisNIDAQ(dirs, DRF.expParam, DRF.DAQin);

dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [AVar.participant AVar.run 'ResultsDRF.mat']);
if AVar.debug == 0
    fprintf('Saving Results for %s %s\n', AVar.participant, AVar.run)
    save(dirs.SavResultsFile, 'auAn', 'auRes', 'niAn', 'niRes')
end
end