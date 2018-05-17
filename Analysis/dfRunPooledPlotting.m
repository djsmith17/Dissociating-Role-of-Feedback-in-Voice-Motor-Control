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
PolPlt.MaskVVoice        = 0;
PolPlt.AllSubjMaskvVoice = 1;

fStat    = 0;

ppi        = 300;
scRes      = [2560 1440];
scDim      = [18.625 11.75];
targFigDim = [15 4];

targPixDim = calcFigPixDim(ppi, scRes, scDim, targFigDim);

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
        drawMeanSubjf0Resp(pooledRunStr(ii), targPixDim, dirs.SavResultsDir, fLabel, fStat, 1)
    end
end

if PolPlt.AllSubjMaskvVoice == 1
    fLabel = 0;
    drawMeanSubjf0Resp(allSubjRes, targPixDim, dirs.SavResultsDir, fLabel, fStat, 1)
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