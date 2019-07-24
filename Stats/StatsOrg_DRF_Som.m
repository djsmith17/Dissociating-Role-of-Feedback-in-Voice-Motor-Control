function StatsOrg_DRF_Som(dirs, pA, allSubjRes)

allSubjStatTable = allSubjRes.statTable;
meas   = {'tAtMin', 'StimMag', 'RespMag', 'RespPer'};
mUnits = {'s', 'cents', 'cents', '%'};

cond    = pA.cond;
numCond = pA.numCond;
pubCond = pA.pubCond;

pubTable = initPubTable(meas, pubCond);
dirs.behavioralResultTable = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'BehavioralResultTable.xlsx']);

curTestingMeas = 1:4;
ApplyTrans = 0;
for k = curTestingMeas
    
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
        measure = curStatTable.(curCond);
        
        measureVar.varName   = meas{k};
        measureVar.condition = curCond;
        measureVar.units     = mUnits{k};
        
        % Perform Standard Summary Stats
        summaryStat = MeasureSummaryStats(dirs, pA, measureVar, measure, lambdas(i));
             
        % Describe the normality
        summaryStat = summaryStat.testNormality();
        
        % Answering Research Questions 1-3
        if k == 3 || k == 4
            summaryStat = summaryStat.performTTest(1); 
            rangeVal = ['A' num2str(7 +2*(i))];
            writetable(summaryStat.statSentTable, dirs.behavioralResultTable, 'Range', rangeVal, 'WriteRowNames', 1, 'Sheet', meas{k})
        end
        
        % Concatenate the Summary Stat Arrays across condition
        summaryVarTableAcrossCond = [summaryVarTableAcrossCond; summaryStat.SummaryTable];

        % Concatenate the Structure for Histogram and Transformed Values
        measureSummaryStrs = cat(1, measureSummaryStrs, summaryStat.SummaryStruct);      
    end
    
    % Find the difference between the two conditions and place in Struct
    measDiff = measureSummaryStrs(1).measure - measureSummaryStrs(2).measure;
    
    measureDiffVar.varName   = [meas{k} 'Diff'];
    measureDiffVar.condition = 'Diff';
    measureDiffVar.units     = mUnits{k};
    summaryStatDiff = MeasureSummaryStats(dirs, pA, measureDiffVar, measDiff, 0);
    
    if k == 4
        % Apply a Box Cox transform with the default (best) lambda
        summaryStatDiff = performSimpleBoxCoxTrans(summaryStatDiff);
    end
    
    summaryStatDiff = summaryStatDiff.testNormality();       % Test Normality
    summaryStatDiff = summaryStatDiff.performTTest();        % Perform t-test
    summaryStatDiff.drawHistoBoxCombo()                      % Visualize Normality/Outliers
    
    % Visualizations
    drawHistograms(measureSummaryStrs, dirs, pA)             % Visualize Distribution/Normality
    drawBoxPlot(measureSummaryStrs, summaryStatDiff.SummaryStruct, dirs, pA)% Visualize Distribution/Outliers 
        
    % Save Behavioral Result Table: Values ready for inclusion in manuscript 
    writetable(summaryVarTableAcrossCond, dirs.behavioralResultTable, 'WriteRowNames', 1, 'Sheet', meas{k})
    writetable(summaryStatDiff.statSentTable, dirs.behavioralResultTable, 'Range', 'A7', 'WriteRowNames', 1, 'Sheet', meas{k})
    
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

function drawHistograms(measureSummaryStrs, dirs, pA)

colors = ['b', 'r', 'g'];
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

dirs.DistributionFigureFile = fullfile(dirs.SavResultsDir, [pAnalysis varName suffix 'DistributionPlot.jpg']);
export_fig(dirs.DistributionFigureFile)
end

function drawBoxPlot(measureSummaryStrs, summaryStrDiff, dirs, pA)

fontN = 'Arial';
axisLSize = 25;

pAnalysis = pA.pAnalysis;
cond      = pA.pubCond;

isSig = summaryStrDiff.isSig;

pAnalysisFix = pAnalysis;
pAnalysisFix(strfind(pAnalysisFix, '_')) = '';

measBox = figure('Color', [1 1 1]);
plotpos = [30 0]; plotdim = [700 1000];
set(measBox, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

varName     = measureSummaryStrs.varName;
measureData = [measureSummaryStrs.measure];

boxplot(measureData, 'Labels', cond)
ylabel([varName ' (' summaryStrDiff.units ')'])
title({pAnalysisFix, varName})
box off

yt = get(gca, 'YTick');
axis([xlim    0  max(yt)*1.2])
xt = get(gca, 'XTick');

if isSig
    hold on
    plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k')
    plot(mean(xt([1 2])), max(yt)*1.12, '*k', 'MarkerSize', 10)
    hold off
end
text(mean(xt([1 2]))-0.25, max(yt)*1.15, ['p = ' summaryStrDiff.ttestPstr], 'FontSize',18)

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold')

dirs.BoxPlotFigureFile = fullfile(dirs.SavResultsDir, [pAnalysis varName 'BoxPlot.jpg']);
export_fig(dirs.BoxPlotFigureFile)
end

function pubTable = initPubTable(meas, pubCond)

numMeas = length(meas);
numCond = length(pubCond);

genVar = cell(numCond, 1);
genVar(:) = {''};

pubTable = table(genVar, genVar, genVar, genVar); %Three times for numMeas
pubTable.Properties.VariableNames = meas;
pubTable.Properties.RowNames = pubCond;
end

function pubTable = popPubTable(pubTable, curCol, summaryVarTableAcrossCond)

[numCond, ~] = size(summaryVarTableAcrossCond);

for ii = 1:numCond
   curMean  = summaryVarTableAcrossCond.mean(ii); % Mean
   curError = summaryVarTableAcrossCond.SD(ii);   % Standard Error of the Mean
   
   curPubPrint = sprintf('%s (%s)', num2str(curMean), num2str(curError));
   pubTable(ii, curCol) = {curPubPrint};
end
end