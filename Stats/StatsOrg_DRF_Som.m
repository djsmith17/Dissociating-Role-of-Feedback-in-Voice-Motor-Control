function StatsOrg_DRF_Som(dirs, pA, allSubjRes)

allSubjStatTable = allSubjRes.statTable;
meas    = {'f0', 'tAtMin', 'StimMag', 'RespMag', 'RespPer'};
measPub = {'Fundamental Frequency', 'Response Latency', 'Stimulus Magnitude', 'Response Magnitude', 'Response Percentage'};
mUnits  = {'Hz', 's', 'cents', 'cents', '%'};
numMeas = length(meas);

measFor1SampleTest  = {'RespMag', 'RespPer'};
measToBeTransformed = {'f0'};

cond    = pA.cond;
numCond = pA.numCond;
pubCond = pA.pubCond;

pubTable = initPubTable(meas, pubCond);
dirs.behavioralResultTable = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'BehavioralResultTable.xlsx']);

ApplyTrans = 1;
for k = 1:numMeas
    
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
        
        % Setup a structure to be fed into the class.
        measureVar.varName    = meas{k};
        measureVar.varNamePub = measPub{k};
        measureVar.condition  = curCond;
        measureVar.units      = mUnits{k};
        
        % Create an object which will make performing stats bvery easy. 
        % Perform Standard Summary Stats
        summaryStat = MeasureSummaryStats(dirs, pA, measureVar, measure, lambdas(i));
             
        % Perform method of class which tests for normality
        summaryStat = summaryStat.testNormality();
        
        % Do we need the result of a 1-sample t-test? 
        % (Signficantly different than 0)
        if ismember(meas{k}, measFor1SampleTest)
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
    
    measureDiffVar.varName     = [meas{k} 'Diff'];
    measureDiffVar.varNamePub  = [meas{k} 'Diff'];
    measureDiffVar.condition   = 'Diff';
    measureDiffVar.units       = mUnits{k};
    summaryStatDiff = MeasureSummaryStats(dirs, pA, measureDiffVar, measDiff, 0);
    
    % Do we need to apply a transformation of the data?
    if ismember(meas{k}, measToBeTransformed) && ApplyTrans
        % Apply a Box Cox transform with the (best) lambda
        summaryStatDiff = summaryStatDiff.performSimpleBoxCoxTrans();
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
fullResultFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'AllVariableTable.csv']);
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

fontN = 'Times New Roman';
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
varNamePub  = measureSummaryStrs.varNamePub;
measureData = [measureSummaryStrs.measure];

boxplot(measureData)
ylabel([varNamePub ' (' summaryStrDiff.units ')'])
% title({pAnalysisFix, varNamePub})
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

if str2num(summaryStrDiff.ttestPstr) < 0.001
    pValSent = 'p < .001';
else
    pValSent = sprintf('p = %0.3f', str2num(summaryStrDiff.ttestPstr));
end

text(mean(xt([1 2]))-0.25, max(yt)*1.15, pValSent, 'FontSize',18)

set(gca, 'XTickLabel', cond)
fix_xticklabels(gca, 0.1, {'FontSize', 17, 'FontName', fontN, 'FontWeight','bold'});

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

pubTable = table(genVar, genVar, genVar, genVar, genVar); % Four times for numMeas
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