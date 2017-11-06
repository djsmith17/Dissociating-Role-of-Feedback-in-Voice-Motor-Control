function dfPooledPlotting()
%The function that draws the plots from analyzed results

clear all; close all; clc
PolPlt.project  = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
PolPlt.poolA    = 'Pooled Analyses'; %List of multiple participants.
PolPlt.analyses = 'SfN2017';

%Plot Toggles. This could eventually become an input variable
sv2File                     = 1;
PolPlt.NIDAQ_MeanTrialMicf0 = 0;
PolPlt.MaskVVoice           = 1;
PolPlt.dataTable            = 0;

dirs                = dfDirs(PolPlt.project);
dirs.SavResultsDir  = fullfile(dirs.Results, PolPlt.poolA, PolPlt.analyses); %Where to save results
dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [PolPlt.analyses 'ResultsDRF.mat']); %Where to save results

if exist(dirs.SavResultsFile, 'file') == 0
    fprintf('\nERROR: File %s does not exist!\n', dirs.SavResultsFile)
    return
end

load(dirs.SavResultsFile)

if PolPlt.NIDAQ_MeanTrialMicf0 == 1
    drawDAQMeanTrialMicf0(combDataStr(1,1), dirs.SavResultsDir)
    drawDAQMeanTrialMicf0(combDataStr(1,2), dirs.SavResultsDir)
    drawDAQMeanTrialMicf0(combDataStr(2,1), dirs.SavResultsDir)
    drawDAQMeanTrialMicf0(combDataStr(2,2), dirs.SavResultsDir)
    drawDAQMeanTrialMicf0(combDataStr(3,1), dirs.SavResultsDir)
    drawDAQMeanTrialMicf0(combDataStr(3,2), dirs.SavResultsDir)
    drawDAQMeanTrialMicf0(combDataStr(4,1), dirs.SavResultsDir)
    drawDAQMeanTrialMicf0(combDataStr(4,2), dirs.SavResultsDir)
end

if PolPlt.MaskVVoice == 1
    drawMaskvVoiceMeanf0(combDataStr(1,1), combDataStr(1,2), statLib(1,:), dirs.SavResultsDir)
    drawMaskvVoiceMeanf0(combDataStr(2,1), combDataStr(2,2), statLib(2,:), dirs.SavResultsDir) 
    drawMaskvVoiceMeanf0(combDataStr(3,1), combDataStr(3,2), statLib(3,:), dirs.SavResultsDir) 
    drawMaskvVoiceMeanf0(combDataStr(4,1), combDataStr(4,2), statLib(4,:), dirs.SavResultsDir) 
end

end