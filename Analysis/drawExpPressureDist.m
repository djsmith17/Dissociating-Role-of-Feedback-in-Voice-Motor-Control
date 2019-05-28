function drawExpPressureDist(dirs, pA, pooledRunStr)

numSubj = length(pooledRunStr);

fontN = 'Arial';
axisLSize = 14;

mu = '\mu';
sigma  = '\sigma'; 

plotpos = [10 10];
plotdim = [1600 800];
IndivTrialAudResp = figure('Color', [1 1 1]);
set(IndivTrialAudResp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

ha = tight_subplot(2, 6, [0.13 0.03],[0.12 0.18],[0.05 0.03]);

allMeanPress = [];
for ii = 1:numSubj
     axes(ha(ii))
     
     curRes  = pooledRunStr(ii);
     curSubj = curRes.studyID;
     curSubj(strfind(curSubj, '_')) = ' ';
     
     pressureOn   = curRes.sensorPOnOff(:,1);
     numTrial     = length(pressureOn);
     pressureOnM  = round(mean(pressureOn), 2);
     pressureOnSD = round(std(pressureOn), 2);
     
     pressureOnSE = pressureOnSD/sqrt(numTrial);
     
     respPerMin = min(pressureOn);
     respPerMax = max(pressureOn);
     lB = respPerMin - 0.02;
     rB = respPerMax + 0.02;
     
     allMeanPress = cat(1, allMeanPress, pressureOnM);
     
     histogram(pressureOn, 10)
     box off
     
     if ii == 7
         xlabel('Pressure (psi)')
     end
     
     set(gca, 'FontName', fontN,...
              'FontSize', axisLSize,...
              'FontWeight','bold')
     
     
     title({[curSubj ' (n=' num2str(numTrial) ')'],[mu '=' num2str(pressureOnM) ', ' sigma '=' num2str(pressureOnSD)]})

     axis([lB rB 0 5])
end

bigTit = suptitle({'Masking Study', 'Distributions of pressure readings at the appex of the pressure step function', [mu '= Mean, ' sigma '= Standard Deviation']});

set(bigTit, 'FontName', fontN,...
              'FontSize', axisLSize,...
              'FontWeight','bold')

pltTitle = 'ExpPressureDistAllSubj.png';
plotFileName = fullfile(dirs.SavResultsDir, pltTitle);
export_fig(plotFileName)

summaryVarTable = RawSummaryStats('Mean Pressure Levels', allMeanPress, 0);
drawHistoBoxCombo(summaryVarTable, dirs, pA)
end

function summaryVarTable = RawSummaryStats(variableName, measureVals, idealLambda)

numObs    = length(measureVals);

% Calculate the Descriptive Stats
summaryVarTable = table();
summaryVarTable.varName = variableName;
summaryVarTable.measure = {measureVals};        % Raw Data Values

summaryVarTable.mean     = round(mean(measureVals), 2);
summaryVarTable.median   = round(median(measureVals), 2);
summaryVarTable.min      = round(min(measureVals), 2);
summaryVarTable.max      = round(max(measureVals), 2);
summaryVarTable.SD       = round(std(measureVals), 2);
summaryVarTable.SE       = round(summaryVarTable.SD/sqrt(numObs), 2);

summaryVarTable.isTrans     = 0;             % Default is not transformed
summaryVarTable.measureT    = {measureVals}; % Transformed Data Values (Default is the same)
summaryVarTable.measureZ    = {measureVals}; % Z-Scored Data Values
summaryVarTable.idealLambda = idealLambda;
summaryVarTable.usedLambda  = 'N/A';   % Default is not transformed
summaryVarTable.suffix      = ' ';      % Default is not transformed

summaryVarTable = testNormality(summaryVarTable);
end

function summaryVarTable = testNormality(summaryVarTable)

% Skew and Kurtosis
summaryVarTable.measureSkew     = round(skewness(summaryVarTable.measureT{1}), 4);
summaryVarTable.measureKurtosis = round(kurtosis(summaryVarTable.measureT{1}), 2);

% Z-Score and Shapiro-Wilk Test
summaryVarTable.measureZ   = {zscore(summaryVarTable.measureT{1})};
[swH, swPValue, swTest]    = swtest(summaryVarTable.measureZ{1});

summaryVarTable.swH      = double(swH);
summaryVarTable.swPValue = round(swPValue, 3);
summaryVarTable.swTest   = round(swTest, 3);
end

function drawHistoBoxCombo(summaryVarTable, dirs, pA)

measure = summaryVarTable.measureT{1};
varName = summaryVarTable.varName;
suffix  = summaryVarTable.suffix;
swH = summaryVarTable.swH; swP = summaryVarTable.swPValue; swW = summaryVarTable.swTest;

pAnalysis = pA.pAnalysis;
mu = '\mu';
sigma  = '\sigma'; 

diffBox = figure('Color', [1 1 1]);
plotpos = [30 0]; plotdim = [800 300];
set(diffBox, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

subplot(1,2,1); histogram(measure, 10); box off
title(['H=' num2str(swH) ', p=' num2str(round(swP,4)) ', W=' num2str(round(swW,3))])

subplot(1,2,2); boxplot(measure); box off
suptitle(varName)

annotation('textbox',[0.80 0.48 0.45 0.1],...
           'string', {[mu ' = ' num2str(summaryVarTable.mean) 'psi'],...
                      [sigma ' = ' num2str(summaryVarTable.SD) 'psi']},...
           'LineStyle','none',...
            'FontWeight','bold',...
            'FontSize',14,...
            'FontName','Arial');

dirs.BoxPlotFigureFile = fullfile(dirs.SavResultsDir, [pAnalysis varName suffix 'BoxPlotCombo.png']);
export_fig(dirs.BoxPlotFigureFile)
end