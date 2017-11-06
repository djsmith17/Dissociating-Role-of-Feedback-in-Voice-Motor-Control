function dfRunPooledAnalysis()
%Analyses the microphone data from the somatosensory perturbation
%experiment. Measures the change in f0 over each trial, and each run for a
%given participant. At the end it approximates a general response to
%inflation to be used in the auditory perturbation experiment

%Require the Signal Processing Toolbox

clear all; close all; clc
pA.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control'; 
pA.pAnalysis     = 'SfN2017';
pA.participants  = {'Pilot24'; 'Pilot25'; 'Pilot26'; 'Pilot22'}; %List of multiple participants.
pA.numPart       = length(pA.participants);
pA.runs          = {'SF1'; 'SF2'; 'SF3'; 'SF4'}; %All runs to consider 
pA.numRuns       = length(pA.runs);

dirs               = dfDirs(pA.project);
dirs.SavResultsDir = fullfile(dirs.Results, 'Pooled Analyses', pA.pAnalysis);

if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

for ii = 1:pA.numPart
    for jj = 1:pA.numRuns
        participant = pA.participants{ii};
        run         = pA.runs{jj};
        dirs.SavFileDir  = fullfile(dirs.Results, participant, run); %Where to save results        
        dirs.SavFile     = fullfile(dirs.SavFileDir, [participant run 'ResultsDRF.mat']);
        
        if exist(dirs.SavFile, 'file') == 0
            disp('ERROR: NO DANG FILE')
            return
        end        
        load(dirs.SavFile)
        AudFB = auAn.AudFB;
        
        if strcmp(AudFb, 'Masking Noise', 2)
           
        else
            
        end
        
        
    end
end

niPAn  = [];
niPRes = [];
auPAn  = [];
auPRes = [];

dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'ResultsDRF.mat']);
fprintf('Saving Pooled Analysis for %s\n', pA.pAnalysis)
save(dirs.SavResultsFile, 'auPAn', 'auPRes', 'niPAn', 'niPRes')
end