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
PolPlt.AllSubjMaskvVoice    = 1;
PolPlt.dataTable            = 0;

dirs                = dfDirs(PolPlt.project);
dirs.SavResultsDir  = fullfile(dirs.Results, PolPlt.poolA, PolPlt.analyses); %Where to save results
dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [PolPlt.analyses 'ResultsDRF.mat']); %Where to save results

ppi        = 300;
scRes      = [1680 1050];
scDim      = [18.625 11.75];
targFigDim = [14 5];

targPixDim = calcFigPixDim(ppi, scRes, scDim, targFigDim);

pltLetter = {'a'; 'b'; 'c'; 'd'};

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
    for ii = 1:4
        pltName = ['SfN2017Results Figure 4' pltLetter{ii}];
        drawMaskvVoiceMeanf0(combDataStr(ii,1), combDataStr(ii,2), statLib(ii,:), targPixDim, pltName, dirs.SavResultsDir)
    end
end

if PolPlt.AllSubjMaskvVoice == 1
    pltName = 'SfN2017Results Figure 5';
    drawMeanSubjf0Resp(allSubjRes, statLibAll, targPixDim, pltName, dirs.SavResultsDir)
end
end

function targPixDim = calcFigPixDim(ppi, scRes, scDim, targFigDim)

mHppi = round(scRes(1)/scDim(1));
mVppi = ceil(scRes(2)/scDim(2));

scaleDif = ppi/mHppi;

tHPixDim = ppi*targFigDim(1);
tVPixDim = ppi*targFigDim(2);

targPixDim = [tHPixDim tVPixDim];
end