function dfRunPooledPlotting()
% dfRunPooledPlotting() loads a pooled results structure and plots these
% pooled results. This function should be used specifically for pooled 
% results and not individual participant results. 
%
% There are a few different plotting functions available:
% -drawMeanTrialMicf0
% -drawMaskvVoiceMeanf0
% -drawMeanSubjf0Resp

close all
PolPlt.project  = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
PolPlt.analyses = 'SfN2017';

dirs                = dfDirs(PolPlt.project);
dirs.SavResultsDir  = fullfile(dirs.Results, 'Pooled Analyses', PolPlt.analyses);       % Analyzed Results Folder
dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [PolPlt.analyses 'ResultsDRF.mat']); % The results to load

% Plot Toggles. Which plots do you want?
PolPlt.MeanTrialMicf0    = 0;
PolPlt.MaskVVoice        = 1;
PolPlt.AllSubjMaskvVoice = 1;

ppi        = 300;
scRes      = [2560 1440];
scDim      = [18.625 11.75];
targFigDim = [15 4];

targPixDim = calcFigPixDim(ppi, scRes, scDim, targFigDim);

pltLetter = {'a'; 'b'; 'c'; 'd'};

if exist(dirs.SavResultsFile, 'file') == 0
    fprintf('\nERROR: File %s does not exist!\n', dirs.SavResultsFile)
    return
else
    load(dirs.SavResultsFile)
    % Returns combDataStr; statLib
end

if PolPlt.MeanTrialMicf0 == 1
    drawMeanTrialMicf0(combDataStr(1,1), dirs.SavResultsDir)
    drawMeanTrialMicf0(combDataStr(1,2), dirs.SavResultsDir)
    drawMeanTrialMicf0(combDataStr(2,1), dirs.SavResultsDir)
    drawMeanTrialMicf0(combDataStr(2,2), dirs.SavResultsDir)
    drawMeanTrialMicf0(combDataStr(3,1), dirs.SavResultsDir)
    drawMeanTrialMicf0(combDataStr(3,2), dirs.SavResultsDir)
    drawMeanTrialMicf0(combDataStr(4,1), dirs.SavResultsDir)
    drawMeanTrialMicf0(combDataStr(4,2), dirs.SavResultsDir)
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
close all
end

function targPixDim = calcFigPixDim(ppi, scRes, scDim, targFigDim)

mHppi = round(scRes(1)/scDim(1));
mVppi = ceil(scRes(2)/scDim(2));

scaleDif = ppi/mHppi;

tHPixDim = ppi*targFigDim(1);
tVPixDim = ppi*targFigDim(2);

targPixDim = [tHPixDim tVPixDim];
end