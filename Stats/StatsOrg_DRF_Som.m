function StatsOrg_DRF_Som(dirs, pA, allSubjRes)

allSubjStatTable = allSubjRes.statTable;
meas = {'StimMag', 'RespMag', 'RespPer'};

cond    = pA.cond;
numCond = pA.numCond;
pubCond = pA.pubCond;

pubTable = initPubTable(meas, pubCond);
dirs.behavioralResultTable = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'BehavioralResultTable.xlsx']);

curTestingMeas = 1:3;
ApplyTrans = 0;
for k = curTestingMeas
    pA.k = k;
    [curStatTable, cond_table] = organizeVarByCond(allSubjStatTable, meas{k}, cond);

    lambdas = [];
    if k == 1 && ApplyTrans
        for i = 1:numCond
            % Identify the Variable and Condition
            measure   = curStatTable.(cond_table{i});

            [~, lambda] = boxcox(measure + 1 - min(measure));
            lambdas = cat(1, lambdas, lambda);
        end
    else
        lambdas = [0 0 0];
    end
     
    measureSummaryStrs        = [];
    summaryVarTableAcrossCond = table();
    for i = 1:numCond
        % Identify the Variable and Condition
        curCond = cond_table{i};
        measure   = curStatTable.(curCond);
        
        % Perform Standard Sumamry Stats
        [summaryVarStr, summaryVarTable] = RawSummaryStats(meas{k}, curCond, measure, lambdas(i));
             
        % Use some function to describe the normality
        summaryVarStr = testNormality(summaryVarStr);
        
        % Add the Normality results to the variableStat Array
        summaryVarTable = concateAdditionalStat(summaryVarStr, summaryVarTable);
        
        % Concatenate the Summary Stat Arrays across condition
        summaryVarTableAcrossCond = [summaryVarTableAcrossCond; summaryVarTable];

        % Concatenate the Structure for Histogram and Transformed Values
        measureSummaryStrs = cat(1, measureSummaryStrs, summaryVarStr);      
    end
    
    % Find the difference between the two conditions and place in Struct
    measDiff = measureSummaryStrs(1).measure - measureSummaryStrs(2).measure;
    [summaryStrDiff, ~] = RawSummaryStats([meas{k} 'Diff'], 'Diff', measDiff, 0);
    
    if k == 3 
        [summaryStrDiff.measureT, l] = boxcox(summaryStrDiff.measure + 1 - min(summaryStrDiff.measure));
        summaryStrDiff.isTrans    = 1;
        summaryStrDiff.suffix     = 'Trans';
        summaryStrDiff.usedLambda = num2str(round(l,2));
    end
    
    % Test for Normality
    summaryStrDiff = testNormality(summaryStrDiff);
    
    % Perform a One-Sample T-Test on the difference between the measures
    [summaryStrDiff.fH, summaryStrDiff.fP] = ttest(summaryStrDiff.measureT);
    summaryStrDiff.fPround = sprintf('%0.6f', summaryStrDiff.fP);
    if summaryStrDiff.fP < (0.05/3)
        summaryStrDiff.isSig = 1;
    else
        summaryStrDiff.isSig = 0;
    end
    
    % Visualizations
    drawHistograms(measureSummaryStrs, dirs, pA)             % Visualize Distribution/Normality
    drawBoxPlot(measureSummaryStrs, summaryStrDiff, dirs, pA)% Visualize Distribution/Outliers
    drawHistoBoxCombo(summaryStrDiff, dirs, pA)              % Visualize Normality/Outliers
        
    % Save Behavioral Result Table: Values ready for inclusion in manuscript 
    writetable(summaryVarTableAcrossCond, dirs.behavioralResultTable, 'WriteRowNames', 1, 'Sheet', meas{k})
    
    % Add to the Table for publication
    pubTable = popPubTable(pubTable, k, summaryVarTableAcrossCond);
end

writetable(pubTable, dirs.behavioralResultTable, 'WriteRowNames', 1, 'Sheet', 'PubTable')

% Save All Variable Table: Easy access excel file to check analysis
fullResultFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'AllVariableTable.xlsx']);
writetable(allSubjStatTable, fullResultFile, 'WriteVariableNames', true)
end

function [curStatTable, cond_Table] = organizeVarByCond(allSubjStatTable, meas, cond)

cond_Table = matlab.lang.makeValidName(cond); % Valid idenitifer strings
condSubj = ['SubjID' cond_Table];

curStatTable = allSubjStatTable(:, {'SubjID', 'AudFB', meas});
curStatTable = unstack(curStatTable, meas, 'AudFB');
curStatTable = curStatTable(:, condSubj);
end

function [summaryVarStr, summaryVarTable] = RawSummaryStats(variableName, cond, measure, idealLambda)

numObs    = length(measure);

% Calculate the Descriptive Stats
summaryVarStr.varName  = variableName;
summaryVarStr.cond     = cond;
summaryVarStr.measure  = measure;        % Raw Data Values

summaryVarStr.mean     = round(mean(measure), 2);
summaryVarStr.median   = round(median(measure), 2);
summaryVarStr.min      = round(min(measure), 2);
summaryVarStr.max      = round(max(measure), 2);
summaryVarStr.SD       = round(std(measure), 2);
summaryVarStr.SE       = round(summaryVarStr.SD/sqrt(numObs), 2);

summaryVarStr.isTrans     = 0;       % Default is not transformed
summaryVarStr.measureT    = measure; % Transformed Data Values (Default is the same)
summaryVarStr.measureZ    = [];      % Z-Scored Data Values
summaryVarStr.idealLambda = idealLambda;
summaryVarStr.usedLambda  = 'N/A';   % Default is not transformed
summaryVarStr.suffix      = '';      % Default is not transformed

summaryVarTable        = table();
summaryVarTable.mean   = summaryVarStr.mean;
summaryVarTable.min    = summaryVarStr.min;
summaryVarTable.median = summaryVarStr.median;
summaryVarTable.max    = summaryVarStr.max;
summaryVarTable.SD     = summaryVarStr.SD;
summaryVarTable.SE     = summaryVarStr.SE;
summaryVarTable.Properties.RowNames = {cond};
end

function summaryStr = testNormality(summaryStr)

% Skew and Kurtosis
summaryStr.measureSkew     = round(skewness(summaryStr.measureT), 4);
summaryStr.measureKurtosis = round(kurtosis(summaryStr.measureT), 2);

% Z-Score and Shapiro-Wilk Test
summaryStr.measureZ        = zscore(summaryStr.measureT);
[swH, swPValue, swTest]    = swtest(summaryStr.measureZ);

summaryStr.swH      = double(swH);
summaryStr.swPValue = round(swPValue, 3);
summaryStr.swTest   = round(swTest, 3);
end

function summaryVarTable = concateAdditionalStat(summaryStr, summaryVarTable)

summaryVarTable.Skew     = summaryStr.measureSkew;
summaryVarTable.Kurtosis = summaryStr.measureKurtosis;
summaryVarTable.swH      = summaryStr.swH;
summaryVarTable.swPValue = summaryStr.swPValue;
summaryVarTable.swTest   = summaryStr.swTest;
end

function drawHistograms(measureSummaryStrs, dirs, pA)

units  = {'cents', 'cents', '%'};
colors = ['b', 'r', 'g'];
sigma  = '\sigma'; mu = '\mu';
lambda = '\lambda';

pAnalysis = pA.pAnalysis;
cond      = pA.cond;
numCond   = pA.numCond;

pAnalysisFix = pAnalysis;
pAnalysisFix(strfind(pAnalysisFix, '_')) = '';

measDist = figure('Color', [1 1 1]);
plotpos = [10 10]; plotdim = [1300 800];
set(measDist, 'Position',[plotpos plotdim],'PaperPositionMode','auto')
ha = tight_subplot(2, numCond,[0.1 0.05],[0.1 0.1],[0.05 0.05]);

for ii = 1:numCond
    summaryStr = measureSummaryStrs(ii);
    
    varName = summaryStr.varName; % Should be the same each time, but this is easier
    suffix  = summaryStr.suffix;  % Should be the same each time, but this is easier
    
    step = 25;
    minBound = floor(summaryStr.min/step)*step;
    maxBound = ceil(summaryStr.max/step)*step;
    distBin  = minBound:step:maxBound;
    nBins    = length(distBin)-1;
    
    if nBins < 6
        nBins = 8;
    end

    axes(ha(ii))
    if summaryStr.isTrans
        histogram(summaryStr.measureT, 8, 'FaceColor', colors(ii), 'EdgeColor', colors(ii))
    else
        histogram(summaryStr.measureT, nBins, 'FaceColor', colors(ii), 'EdgeColor', colors(ii), 'BinLimits', [minBound maxBound])
    end
    
    histoTitle = [cond{ii} ' (WStat = ' num2str(summaryStr.swTest) ')'];
    title(histoTitle)
    box off;

    axes(ha(ii+numCond))
    cdfplot(summaryStr.measureZ)
    hold on
    xValues = linspace(min(summaryStr.measureZ), max(summaryStr.measureZ));
    plot(xValues, normcdf(xValues, 0, 1), 'r-')
    legend('Empirical CDF','Standard Normal CDF','Location','best') 
end
suptitle({pAnalysisFix, [varName suffix]})

annotation('textbox',[0.8 0.88 0.45 0.1],...
           'string', {[lambda ' = ' summaryStr.usedLambda]},...
           'LineStyle','none',...
            'FontWeight','bold',...
            'FontSize',14,...
            'FontName','Arial');

dirs.DistributionFigureFile = fullfile(dirs.SavResultsDir, [pAnalysis varName suffix 'DistributionPlot.png']);
export_fig(dirs.DistributionFigureFile)
end

function drawBoxPlot(measureSummaryStrs, summaryStrDiff, dirs, pA)

units  = {'cents', 'cents', '%'};
fontN = 'Arial';
axisLSize = 25;

pAnalysis = pA.pAnalysis;
cond      = pA.pubCond;
numCond   = pA.numCond;

isSig = summaryStrDiff.isSig;

pAnalysisFix = pAnalysis;
pAnalysisFix(strfind(pAnalysisFix, '_')) = '';

measBox = figure('Color', [1 1 1]);
plotpos = [30 0]; plotdim = [700 1000];
set(measBox, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

collData = [];
for i = 1:numCond
    collData(:,i) = measureSummaryStrs(i).measure;
    varName = measureSummaryStrs(i).varName;
end

boxplot(collData, 'Labels', cond)
ylabel([varName ' (' units{pA.k} ')'])
title({pAnalysisFix, varName})
box off

yt = get(gca, 'YTick');
axis([xlim    0  ceil(max(yt)*1.2)])
xt = get(gca, 'XTick');

if isSig
    hold on
    plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k')
    plot(mean(xt([1 2])), max(yt)*1.12, '*k', 'MarkerSize', 10)
    hold off
end
text(mean(xt([1 2]))-0.25, max(yt)*1.15, ['p = ' summaryStrDiff.fPround], 'FontSize',18)

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold')

dirs.BoxPlotFigureFile = fullfile(dirs.SavResultsDir, [pAnalysis varName 'BoxPlot.png']);
export_fig(dirs.BoxPlotFigureFile)
end

function drawHistoBoxCombo(summaryStr, dirs, pA)

measure = summaryStr.measureT;
varName = summaryStr.varName;
suffix  = summaryStr.suffix;
swH = summaryStr.swH; swP = summaryStr.swPValue; swW = summaryStr.swTest;

pAnalysis = pA.pAnalysis;
lambda = '\lambda';

diffBox = figure('Color', [1 1 1]);
plotpos = [30 0]; plotdim = [800 300];
set(diffBox, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

subplot(1,2,1); histogram(measure, 10); box off
title(['H=' num2str(swH) ', p=' num2str(round(swP,4)) ', W=' num2str(round(swW,3))])

subplot(1,2,2); boxplot(measure); box off
suptitle(varName)

annotation('textbox',[0.8 0.88 0.45 0.1],...
           'string', {[lambda ' = ' summaryStr.usedLambda]},...
           'LineStyle','none',...
            'FontWeight','bold',...
            'FontSize',14,...
            'FontName','Arial');

dirs.BoxPlotFigureFile = fullfile(dirs.SavResultsDir, [pAnalysis varName suffix 'BoxPlotCombo.jpg']);
export_fig(dirs.BoxPlotFigureFile)
end

function pubTable = initPubTable(meas, pubCond)

numMeas = length(meas);
numVar  = length(pubCond);

genVar = {''; ''};

pubTable = table(genVar, genVar, genVar);
pubTable.Properties.VariableNames = meas;
pubTable.Properties.RowNames = pubCond;
end

function pubTable = popPubTable(pubTable, curCol, summaryVarTableAcrossCond)

[numCond, ~] = size(summaryVarTableAcrossCond);

for ii = 1:numCond
   curMean = summaryVarTableAcrossCond.mean(ii);
   curSE   = summaryVarTableAcrossCond.SE(ii);
   
   curPubPrint = sprintf('%s (%s)', num2str(curMean), num2str(curSE));
   pubTable(ii, curCol) = {curPubPrint};
end
end