function dfRunSubjPlotting()
%The function that draws the plots from analyzed results

%This uses the following plot functions
%drawDAQAll
%drawDAQPresMic
%drawDAQAlignedPressure
%drawDAQMeanTrialMicf0
%drawAudRespIndivTrial
%drawAudRespMeanTrial

clear all; close all; clc
sPlt.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
sPlt.participants  = {'PureTone200'}; %List of multiple participants.
sPlt.numPart       = length(sPlt.participants);
sPlt.runs          = {'AF1'}; %All runs to consider 
sPlt.numRuns       = length(sPlt.runs);
dirs               = dfDirs(sPlt.project);

%Plot Toggles. This could eventually become an input variable
sv2File                      = 1;
sPlt.drawDAQAll              = 0; % All signals recorded by the NIDAQ
sPlt.drawDAQPresMic          = 0; % Pressure vs Microphone Data
sPlt.drawDAQAlignedPressure  = 0; % Superimposed Pressure recordings from perturbed trials
sPlt.drawDAQMeanTrialMicf0   = 0; % Mean Trials Microphone input. Control vs Perturbed Trials
sPlt.drawDAQMeanTrialAudResp = 1; % Mean Perturbed Trials. Microphone vs Headphones

sPlt.NIDAQ_AllPertTrial      = 0; % 

for ii = 1:sPlt.numPart
    participant = sPlt.participants{ii};
    for jj = 1:sPlt.numRuns 
        run = sPlt.runs{jj};
        dirs.SavResultsDir  = fullfile(dirs.Results, participant, run); %The Analyzed Results Folder...Where Plots will go
        dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [participant run 'ResultsDRF.mat']); %The Analyzed Results FIle

        if exist(dirs.SavResultsFile, 'file') == 0
            fprintf('\nERROR: File %s does not exist!\n', dirs.SavResultsFile)
            return
        else
            load(dirs.SavResultsFile)
        end        
        
        if sPlt.drawDAQAll == 1
            drawDAQAll(niAn, dirs.SavResultsDir, sv2File)
        end
        
        if sPlt.drawDAQPresMic == 1
            drawDAQPresMic(niRes, dirs.SavResultsDir, sv2File)
        end
        
        if sPlt.drawDAQAlignedPressure == 1
            drawDAQAlignedPressure(niRes, dirs.SavResultsDir, sv2File)
        end
        
        if sPlt.drawDAQMeanTrialMicf0 == 1
            drawDAQMeanTrialMicf0(niRes, dirs.SavResultsDir)
        end
        
        if sPlt.drawDAQMeanTrialAudResp == 1
            drawAudRespIndivTrial(niRes, dirs.SavResultsDir)
            drawAudRespMeanTrial(niRes, dirs.SavResultsDir)
        end
        
              
        if sPlt.drawDAQ == 1
            drawDAQAllPertTrialMicf0(niRes, dirs.SavResultsDir)
        end 
    end
end
end