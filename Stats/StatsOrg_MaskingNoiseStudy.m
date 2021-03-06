function StatsOrg_MaskingNoiseStudy(dirs, pA, allSubjRes)

allSubjStatTable = allSubjRes.statTable;
meas = {'StimMag', 'RespMag', 'RespPer', 'tAtMin'};
mUnits = {'cents', 'cents', '%', 's'};

cond    = pA.cond;
numCond = pA.numCond;
pubCond = pA.pubCond;
alphaLevel    = 0.05/3;

pubTable = initPubTable(meas, pubCond);
dirs.behavioralResultTable = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'BehavioralResultTable.xlsx']);

curTestingMeas = 1:4;
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
        
        measureVar.varName   = meas{k};
        measureVar.condition = curCond;
        measureVar.units     = mUnits{k};
        
        % Perform Standard Summary Stats
        summaryStat = MeasureSummaryStats(dirs, pA, measureVar, measure, lambdas(i));
        
        if k == 1 && ApplyTrans
            usedLambda = lambdas(1);
            summaryStat.SummaryStruct.measureT   = boxcox(usedLambda, summaryStat.SummaryStruct.measure);
            summaryStat.SummaryStruct.isTrans    = 1;
            summaryStat.SummaryStruct.suffix     = 'TransVF';
            summaryStat.SummaryStruct.usedLambda = num2str(usedLambda);
        elseif k == 3 && ApplyTrans
            usedLambda = lambdas(3);
            summaryStat.SummaryStruct.measureT   = boxcox(usedLambda, [summaryStat.SummaryStruct.measure + 1 - summaryStat.SummaryStruct.min]);
            summaryStat.SummaryStruct.isTrans    = 1;
            summaryStat.SummaryStruct.suffix     = 'TransACBC';
            summaryStat.SummaryStruct.usedLambda = num2str(usedLambda);
        end
             
        % Describe the normality
        summaryStat = summaryStat.testNormality();
        
        % Concatenate the Summary Stat Arrays across condition
        summaryVarTableAcrossCond = [summaryVarTableAcrossCond; summaryStat.SummaryTable];

        % Concatenate the Structure for Histogram and Transformed Values
        measureSummaryStrs = cat(1, measureSummaryStrs, summaryStat.SummaryStruct);
    end
    
    % Visualizations
    drawHistograms(measureSummaryStrs, dirs, pA)
    drawBoxPlot(measureSummaryStrs, dirs, pA)
    
    if k == 2
        statSentTable = testParametric(curStatTable, cond_table, meas{k}, alphaLevel);
    else
        statSentTable = testNonParametric(curStatTable, meas{k}, alphaLevel);
    end
    
    % Save Behavioral Result Table: Values ready for inclusion in manuscript
    writetable(summaryVarTableAcrossCond, dirs.behavioralResultTable, 'WriteRowNames', 1, 'Sheet', meas{k})
    writetable(statSentTable, dirs.behavioralResultTable, 'Range', 'A7', 'WriteRowNames', 1, 'Sheet', meas{k})
    
    % Add to the Table for publication
    pubTable = popPubTable(pubTable, k, summaryVarTableAcrossCond);
end

writetable(pubTable, dirs.behavioralResultTable, 'WriteRowNames', 1, 'Sheet', 'PubTable')

% Save All Variable Table: Easy access excel file to check analysis
fullResultFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'AllVariableTable.xlsx']);
writetable(allSubjStatTable, fullResultFile, 'WriteVariableNames',true)
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
cond      = pA.pubCond;
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

dirs.DistributionFigureFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis varName suffix 'DistributionPlot.jpg']);
export_fig(dirs.DistributionFigureFile)
end

function drawBoxPlot(measureSummaryStrs, dirs, pA)

fontN = 'Times New Roman';
axisLSize = 25;

pAnalysis = pA.pAnalysis;
cond      = pA.pubCond;

measBox = figure('Color', [1 1 1]);
plotpos = [30 30]; plotdim = [700 1000];
set(measBox, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

varName     = measureSummaryStrs.varName;
units       = measureSummaryStrs.units;  
measureData = [measureSummaryStrs.measure];

boxplot(measureData)
ylabel([varName ' (' units ')'])
title(varName)
box off

set(gca, 'XTickLabel', cond)
fix_xticklabels(gca, 0.1, {'FontSize', 17, 'FontName', fontN, 'FontWeight','bold'});

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold')

dirs.BoxPlotFigureFile = fullfile(dirs.SavResultsDir, [pAnalysis varName 'BoxPlot.jpg']);
export_fig(dirs.BoxPlotFigureFile)
end

function statSentTable = testParametric(curStatTable, cond_table, varName, alpha)

condTable = table(cond_table');

measFit = fitrm(curStatTable, 'VoiceFeedback-AC_BCMaskingNoise~1', 'WithinDesign', condTable, 'WithinModel', 'separatemeans');
measSph = mauchly(measFit);
    
rAnovaRes = ranova(measFit);
pVal = rAnovaRes.pValue(1);

if pVal < alpha
    sigNote = '';
else
    sigNote = ' not';
end

if measSph.pValue < alpha
    sphNote = ' not';
else
    sphNote = '';
end

statSentence = sprintf('Statistical analyses revealed that there was%s a significant effect of auditory feedback condition on %s (F(%d,%d) = %0.2f, p = %0.3f)\n',...
                       sigNote,...
                       varName,...
                       rAnovaRes.DF(1),...
                       rAnovaRes.DF(2),...
                       rAnovaRes.F(1),...
                       pVal);

spherSentence = sprintf('%s did%s met the assumption of sphericity (chisq(%d) = %0.4f, p = %0.2f)\n', varName, sphNote, measSph.DF, measSph.ChiStat, measSph.pValue);

statSentTable = table({statSentence}, {spherSentence}, 'VariableNames', {'StatSentence', 'SphericitySentence'});
end

function statSentTable = testNonParametric(curStatTable, varName, alpha)

matVer = curStatTable{:,2:4};

[pFried, Ftable, ~] = friedman(matVer, 1, 'off');

if pFried < alpha
    sigNote = '';
else
    sigNote = ' not';
end

statSentence = sprintf('Statistical analyses revealed that there was%s a significant effect of auditory feedback condition on %s (chiSq(%d) = %0.2f, p = %0.3f)\n',...
                       sigNote,...
                       varName,...
                       Ftable{2,3},...
                       Ftable{2,2},...
                       pFried);

statSentTable = table({statSentence}, 'VariableNames', {'StatSentence'});
end

function pubTable = initPubTable(meas, pubCond)

numMeas = length(meas);
numCond  = length(pubCond);

genVar = cell(numCond, 1);
genVar(:) = {''};

pubTable = table(genVar, genVar, genVar, genVar); % Three times for numMeas
pubTable.Properties.VariableNames = meas;
pubTable.Properties.RowNames = pubCond;
end

function pubTable = popPubTable(pubTable, curCol, summaryVarTableAcrossCond)

[numCond, ~] = size(summaryVarTableAcrossCond);

for ii = 1:numCond
   curMean  = summaryVarTableAcrossCond.mean(ii); % Mean
   curError = summaryVarTableAcrossCond.SD(ii);   % Standard Definition
   
   curPubPrint = sprintf('%s (%s)', num2str(curMean), num2str(curError));
   pubTable(ii, curCol) = {curPubPrint};
end
end