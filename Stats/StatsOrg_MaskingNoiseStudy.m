function StatsOrg_MaskingNoiseStudy(dirs, pA, allSubjRes)

allSubjStatTable = allSubjRes.statTable;
meas = {'StimMag', 'RespMag', 'RespPer'};

cond    = pA.cond;
numCond = pA.numCond;

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
     
    measureSummaryStrs     = [];
    variableStatAcrossCond = [];
    for i = 1:numCond
        % Identify the Variable and Condition
        measure   = curStatTable.(cond_table{i});
        
        % Perform Standard Sumamry Stats
        [summaryStr, variableStat] = RawSummaryStats(meas{k}, measure);
        
        summaryStr.idealLambda = lambdas(i);
        
        if k == 1 && ApplyTrans
            usedLambda = lambdas(1);
            summaryStr.measureT   = boxcox(usedLambda, summaryStr.measure);
            summaryStr.isTrans    = 1;
            summaryStr.suffix     = 'TransVF';
            summaryStr.usedLambda = num2str(usedLambda);
        elseif k == 3 && ApplyTrans
            usedLambda = lambdas(3);
            summaryStr.measureT   = boxcox(usedLambda, [summaryStr.measure + 1 - min(summaryStr.measure)]);
            summaryStr.isTrans    = 1;
            summaryStr.suffix     = 'TransACBC';
            summaryStr.usedLambda = num2str(usedLambda);
        else
            summaryStr.measureT   = summaryStr.measure;
            summaryStr.isTrans    = 0;
            summaryStr.suffix     = '   noTrans';
            summaryStr.usedLambda = 'N/A';
        end
             
        % Use some function to describe the normality
        summaryStr = testNormality(summaryStr);
        
        % Add the Normality results to the variableStat Array
        variableStat = concateAdditionalStat(summaryStr, variableStat);
        
        % Concatenate the Summary Stat Arrays across condition
        variableStatAcrossCond = cat(2, variableStatAcrossCond, variableStat);

        % Concatenate the Structure for Histogram and Transformed Values
        measureSummaryStrs = cat(1, measureSummaryStrs, summaryStr);      
    end
    plotHistograms(measureSummaryStrs, dirs, pA)
    drawBoxPlot(measureSummaryStrs, dirs, pA)
    
    if k == 2
        [rAnovaRes, measSph] = testParametric(curStatTable, cond_table);
    else
        [tFried] = testNonParametric(curStatTable);
    end

    dirs.behavioralResultTable = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'BehavioralResultTable' summaryStr.suffix '.xlsx']);
    xlswrite(dirs.behavioralResultTable, variableStatAcrossCond, meas{k})
end

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

function [summaryStr, variableStat] = RawSummaryStats(variableName, measure)

numObs    = length(measure);

% Calculate the Descriptive Stats
summaryStr.varName  = variableName;
summaryStr.measure  = measure;        % Raw Data Values
summaryStr.measureT = [];             % Transformed Data Values
summaryStr.measureZ = [];             % Z-Scored Data Values
summaryStr.idealLambda = [];
summaryStr.usedLambda  = [];
summaryStr.mean     = round(mean(measure), 2);
summaryStr.median   = round(median(measure), 2);
summaryStr.min      = round(min(measure), 2);
summaryStr.max      = round(max(measure), 2);
summaryStr.SD       = round(std(measure), 2);
summaryStr.SE       = round(summaryStr.SD/sqrt(numObs), 2);

variableStat = {summaryStr.mean;...
                summaryStr.median;...
                summaryStr.min;...
                summaryStr.max;...
                summaryStr.SD;...
                summaryStr.SE};
                  

end

function drawBoxPlot(measureSummaryStrs, dirs, pA)

units  = {'cents', 'cents', '%'};
colors = ['b', 'r', 'g'];

cond    = pA.cond;
numCond = pA.numCond;

measBox = figure('Color', [1 1 1]);
plotpos = [30 30]; plotdim = [700 1000];
set(measBox, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

collData = [];
for i = 1:numCond
    collData(:,i) = measureSummaryStrs(i).measure;
    varName = measureSummaryStrs(i).varName;
end

boxplot(collData, 'Labels', cond, 'Whisker', 30)
xlabel('AudFB')
ylabel([varName ' (' units{pA.k} ')'])
title(varName)
box off

dirs.BoxPlotFigureFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis varName 'BoxPlot.jpg']);
export_fig(dirs.BoxPlotFigureFile)
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

function variableStat = concateAdditionalStat(summaryStr, variableStat)

variableStat = cat(1, variableStat, summaryStr.measureSkew);
variableStat = cat(1, variableStat, summaryStr.measureKurtosis);
variableStat = cat(1, variableStat, summaryStr.swH);
variableStat = cat(1, variableStat, summaryStr.swPValue);
variableStat = cat(1, variableStat, summaryStr.swTest);
end

function plotHistograms(measureSummaryStrs, dirs, pA)

units  = {'cents', 'cents', '%'};
colors = ['b', 'r', 'g'];
sigma  = '\sigma'; mu = '\mu';
lambda = '\lambda';

cond    = pA.cond;
numCond = pA.numCond;

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
suptitle({pA.pAnalysis, [varName suffix]})

annotation('textbox',[0.8 0.88 0.45 0.1],...
           'string', {[lambda ' = ' summaryStr.usedLambda]},...
           'LineStyle','none',...
            'FontWeight','bold',...
            'FontSize',14,...
            'FontName','Arial');

dirs.DistributionFigureFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis varName suffix 'DistributionPlot.jpg']);
export_fig(dirs.DistributionFigureFile)
end

function [rAnovaRes, measSph] = testParametric(curStatTable, cond_table)

condTable = table(cond_table');

measFit = fitrm(curStatTable, 'VoiceFeedback-AC_BCMaskingNoise~1', 'WithinDesign', condTable, 'WithinModel', 'separatemeans');
measSph = mauchly(measFit);
    
rAnovaRes = ranova(measFit);
end

function [tFried] = testNonParametric(curStatTable)

matVer = curStatTable{:,2:4};

[pFried,tFried, stats] = friedman(matVer);

end