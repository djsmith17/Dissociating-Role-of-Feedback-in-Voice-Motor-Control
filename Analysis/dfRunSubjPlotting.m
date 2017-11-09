function dfRunSubjPlotting()
%The function that draws the plots from analyzed results

clear all; close all; clc
sPlt.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
sPlt.participants  = {'Pilot24'; 'Pilot25'; 'Pilot26'; 'Pilot22'}; %List of multiple participants.
sPlt.numPart       = length(sPlt.participants);
sPlt.runs          = {'SF1'; 'SF2'; 'SF3'; 'SF4'}; %All runs to consider 
sPlt.numRuns       = length(sPlt.runs);
dirs               = dfDirs(sPlt.project);

%Plot Toggles. This could eventually become an input variable
sv2File              = 1;
sPlt.NIDAQ_allCh     = 0; %Voltage trace of force sensor signal
sPlt.NIDAQ_PresMic   = 0;
sPlt.NIDAQ_AligSens  = 0;
sPlt.NIDAQ_AllPertTrial = 1;
sPlt.NIDAQ_MeanTrialMicf0 = 0;
sPlt.IntraTrial_T    = 0; %SPL trace of individual trial
sPlt.IntraTrial_f0   = 0; %f0 trace for each individual trial
sPlt.InterTrial_f0   = 0; %Average f0 trace over all trials of a run
sPlt.IntraTrialP_f0  = 0; %f0 trace of pertrubed trials of a run
sPlt.InterRun_f0       = 0; %Average f0 trace over all runs analyzed
sPlt.InterTrial_AudRes = 0; %Average f0 response trace to auditory pert trials of a run
sPlt.InterRun_AudRes   = 0; %Average f0 response trace to auditory pert over all runs analyzed

for ii = 1:sPlt.numPart
    participant = sPlt.participants{ii};
    for jj = 1:sPlt.numRuns 
        run = sPlt.runs{jj};
        dirs.SavResultsDir  = fullfile(dirs.Results, participant, run); %Where to save results
        dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [participant run 'ResultsDRF.mat']); %Where to save results

        if exist(dirs.SavResultsFile, 'file') == 0
            fprintf('\nERROR: File %s does not exist!\n', dirs.SavResultsFile)
            return
        end

        %Load 'auAn' 'auRes' 'niAn' 'niRes'
        load(dirs.SavResultsFile)

        if sPlt.NIDAQ_AllPertTrial == 1
            drawDAQAllPertTrialMicf0(niRes, dirs.SavResultsDir)
        end  

        if sPlt.NIDAQ_MeanTrialMicf0 == 1
            drawDAQMeanTrialMicf0(niRes, dirs.SavResultsDir)
        end
        
        if sPlt.NIDAQ_AligSens == 1
            drawDAQAlignedPressure(niRes, dirs.SavResultsDir, sv2File)
        end

        if sPlt.InterTrial_f0 == 1
            drawInterTrialf0(auRes.timeSec, auRes.meanTrialf0_St, auRes.meanTrialf0_Sp, auRes.f0LimitsSec, auRes.trialCount, auRes.meanTrialf0b, auAn.curSess, '', dirs.SavResultsDir)
        end

        if sPlt.NIDAQ_PresMic == 1
            drawDAQPresMic(niRes, dirs.SavResultsDir, sv2File)
        end

        if sPlt.NIDAQ_AligSens == 1
            drawDAQAlignedPressure(niRes, dirs.SavResultsDir, sv2File)
        end

        if sPlt.NIDAQ_allCh == 1
            drawDAQAll(niAn, 2, dirs.SavResultsDir, sv2File)
        end

        if sPlt.IntraTrialP_f0 == 1
            drawAllTrialf0(auRes.time, auRes.allTrialf0, auRes.runTrialOrder, auAn.trigsT, auRes.f0Limits, auRes.meanTrialf0b, auAn.curSess, '', dirs.SavResultsDir)
        end

        if sPlt.InterTrial_AudRes == 1
            drawAudResp_AllTrial(auRes, auAn.curSess, DRF.expParam.curSess, dirs.SavResultsDir)
        end

        if sPlt.InterRun_AudRes == 1
            drawAudResp_InterTrial(auRes.timeSec, auRes.meanTrialf0_St(:,:,2), auRes.meanTrialf0_Sp(:,:,2), auRes.f0LimitsSec, auRes.trialCount, auRes.meanTrialf0b, auAn.curSess, '', dirs.SavResultsDir)
        end  
    end
end
end