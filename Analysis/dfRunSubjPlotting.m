function dfRunSubjPlotting()
%The function that draws the plots from analyzed results

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
sPlt.InterTrial_f0   = 0; %Average f0 trace over all trials of a run
sPlt.IntraTrialP_f0  = 0; %f0 trace of pertrubed trials of a run

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
        

        if sPlt.InterTrial_f0 == 1
            drawInterTrialf0(auRes.timeSec, auRes.meanTrialf0_St, auRes.meanTrialf0_Sp, auRes.f0LimitsSec, auRes.trialCount, auRes.meanTrialf0b, auAn.curSess, '', dirs.SavResultsDir)
        end
        

        if sPlt.IntraTrialP_f0 == 1
            drawAllTrialf0(auRes.time, auRes.allTrialf0, auRes.runTrialOrder, auAn.trigsT, auRes.f0Limits, auRes.meanTrialf0b, auAn.curSess, '', dirs.SavResultsDir)
        end

    end
end
end