function dfRunSubjPlotting()
%The function that draws the plots from analyzed results

clear all; close all; clc
PltVar.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
PltVar.participant  = 'Pilot26'; %List of multiple participants.
PltVar.run          = 'SF4';

%Plot Toggles. This could eventually become an input variable
sv2File                = 0;
PltVar.NIDAQ_allCh     = 0; %Voltage trace of force sensor signal
PltVar.NIDAQ_PresMic   = 0;
PltVar.NIDAQ_AligSens  = 0;
PltVar.NIDAQ_AllPertTrial = 1;
PltVar.NIDAQ_MeanTrialMicf0 = 1;
PltVar.IntraTrial_T    = 0; %SPL trace of individual trial
PltVar.IntraTrial_f0   = 0; %f0 trace for each individual trial
PltVar.InterTrial_f0   = 0; %Average f0 trace over all trials of a run
PltVar.IntraTrialP_f0  = 0; %f0 trace of pertrubed trials of a run
PltVar.InterRun_f0       = 0; %Average f0 trace over all runs analyzed
PltVar.InterTrial_AudRes = 0; %Average f0 response trace to auditory pert trials of a run
PltVar.InterRun_AudRes   = 0; %Average f0 response trace to auditory pert over all runs analyzed

dirs                = dfDirs(PltVar.project);
dirs.SavResultsDir  = fullfile(dirs.Results, PltVar.participant, PltVar.run); %Where to save results
dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [PltVar.participant PltVar.run 'ResultsDRF.mat']); %Where to save results

if exist(dirs.SavResultsFile, 'file') == 0
    fprintf('\nERROR: File %s does not exist!\n', dirs.SavResultsFile)
    return
end

%Load 'auAn' 'auRes' 'niAn' 'niRes'
load(dirs.SavResultsFile)

if PltVar.NIDAQ_AllPertTrial == 1
    drawDAQAllPertTrialMicf0(niRes, dirs.SavResultsDir)
end  

if PltVar.NIDAQ_MeanTrialMicf0 == 1
    drawDAQMeanTrialMicf0(niRes, dirs.SavResultsDir)
end

if PltVar.InterTrial_f0 == 1
    drawInterTrialf0(auRes.timeSec, auRes.meanTrialf0_St, auRes.meanTrialf0_Sp, auRes.f0LimitsSec, auRes.trialCount, auRes.meanTrialf0b, auAn.curSess, '', dirs.SavResultsDir)
end

if PltVar.NIDAQ_PresMic == 1
    drawDAQPresMic(niRes, dirs.SavResultsDir, sv2File)
end

if PltVar.NIDAQ_AligSens == 1
    drawDAQAlignedPressure(niRes, dirs.SavResultsDir, sv2File)
end

if PltVar.NIDAQ_allCh == 1
    drawDAQAll(niAn, 2, dirs.SavResultsDir, sv2File)
end

if PltVar.IntraTrialP_f0 == 1
    drawAllTrialf0(auRes.time, auRes.allTrialf0, auRes.runTrialOrder, auAn.trigsT, auRes.f0Limits, auRes.meanTrialf0b, auAn.curSess, '', dirs.SavResultsDir)
end

if PltVar.InterTrial_AudRes == 1
    drawAudResp_AllTrial(auRes, auAn.curSess, DRF.expParam.curSess, dirs.SavResultsDir)
end

if PltVar.InterRun_AudRes == 1
    drawAudResp_InterTrial(auRes.timeSec, auRes.meanTrialf0_St(:,:,2), auRes.meanTrialf0_Sp(:,:,2), auRes.f0LimitsSec, auRes.trialCount, auRes.meanTrialf0b, auAn.curSess, '', dirs.SavResultsDir)
end         
end