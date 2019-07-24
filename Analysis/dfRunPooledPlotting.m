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
PolPlt.analyses = 'DRF_Som';

dirs                = dfDirs(PolPlt.project);
dirs.SavResultsDir  = fullfile(dirs.Results, 'Pooled Analyses', PolPlt.analyses);       % Analyzed Results Folder
dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [PolPlt.analyses 'ResultsDRF.mat']); % The results to load

% Plot Toggles. Which plots do you want?
PolPlt.MeanTrialMicf0    = 0;
PolPlt.MaskVVoice        = 0;
PolPlt.AllSubjMaskvVoice = 1;
PolPlt.MicVHead          = 0;
PolPlt.AllSubjMicVHead   = 0;

fStat    = 0;
fPres    = 0;

targPPI    = 300;
scRes      = [2560 1440];
scDim      = [6 4]; %inches
monSize    = [23.5 13.5]; %inches ~100ppi
monitorPPI = 100;

targPixDim = calcFigPixDim(targPPI, scDim, monitorPPI);

if exist(dirs.SavResultsFile, 'file') == 0
    fprintf('\nERROR: File %s does not exist!\n', dirs.SavResultsFile)
    return
else
    load(dirs.SavResultsFile)
    % Returns pooledRunStr, allSubjRes
end

if PolPlt.MeanTrialMicf0 == 1
    numIndivi = length(pltNm.pltNameMVi);
    
    for ii = 1:numIndivi
        drawMeanTrialMicf0(combDataStr(ii,1), dirs.SavResultsDir)
        drawMeanTrialMicf0(combDataStr(ii,2), dirs.SavResultsDir)
    end
end

if PolPlt.MaskVVoice == 1
    fLabel = 0;
    numIndivi = length(pooledRunStr);
    for ii = 1:numIndivi
        drawMeanSubjf0Resp(pooledRunStr(ii), targPixDim, dirs.SavResultsDir, fLabel, fStat, fPres)
    end
end

if PolPlt.AllSubjMaskvVoice == 1
    fLabel = 0;
    drawMeanSubjf0Resp_Onset(allSubjRes, targPixDim, dirs.SavResultsDir, fLabel, fStat, fPres)
end

if PolPlt.MicVHead == 1
    fLabel = 0;
    numIndivi = length(pooledRunStr);
    for ii = 1:numIndivi
        drawMeanSubjMicHeadResp(pooledRunStr(ii), targPixDim, dirs.SavResultsDir, fLabel)
    end
end

if PolPlt.AllSubjMicVHead == 1
    fLabel = 0;
    drawMeanSubjMicHeadResp(allSubjRes, targPixDim, dirs.SavResultsDir, fLabel)
end

close all
end

function targPixDim = calcFigPixDim(targPPI, scDim, monitorPPI)

scalePPI = targPPI/monitorPPI;

tHPixDim = scalePPI*(monitorPPI*scDim(1));
tVPixDim = scalePPI*(monitorPPI*scDim(2));

targPixDim = [tHPixDim tVPixDim];
end