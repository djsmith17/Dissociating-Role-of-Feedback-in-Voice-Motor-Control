function StatsOrg_MaskingNoiseStudy(dirs, pA, allSubjRes)

dirs.behavioralResultTable = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'BehavioralResultTable.xlsx']);

allSubjStatTable = allSubjRes.statTable;
meas = {'StimMag', 'RespMag', 'RespPer'};

cond    = pA.cond;
numCond = pA.numCond;

curTestingMeas = 1;
for k = curTestingMeas
    curStatTable = allSubjStatTable(:, {'AudFB', meas{k}});

    measureSummaryStrs     = [];
    measureTAcrossCond     = [];
    variableStatAcrossCond = [];
    for i = 1:numCond
        curFB = strcmp(curStatTable.AudFB, cond(i));
        
        % Identify the Variable and Condition
        measure   = curStatTable{curFB, 2};
        
        % Perform Standard Sumamry Stats
        [summaryStr, variableStat] = RawSummaryStats(meas{k}, measure);
        
        if k == 4
            summaryStr.measureT = boxcox(summaryStr.measure);
            summaryStr.isTrans = 1;
            summaryStr.suffix  = 'Trans';
        else
            summaryStr.measureT = summaryStr.measure;
            summaryStr.isTrans = 0;
            summaryStr.suffix  = '';
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
%     measFit = fitrm(curStatTable, 'RespMag~AudFB');
%     measSph = mauchly(measFit);
    
    xlswrite(dirs.behavioralResultTable, variableStatAcrossCond, meas{k})
end

fullResultFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'AllVariableTable.xlsx']);
writetable(allSubjStatTable, fullResultFile, 'WriteVariableNames',true)
end

function [summaryStr, variableStat] = RawSummaryStats(variableName, measure)

numObs    = length(measure);

% Calculate the Descriptive Stats
summaryStr.varName  = variableName;
summaryStr.measure  = measure;        % Raw Data Values
summaryStr.measureT = [];             % Transformed Data Values
summaryStr.measureZ = [];             % Z-Scored Data Values 
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
mu     = '\mu';
sigma  = '\sigma';

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
    title(cond{ii})
    box off;

    axes(ha(ii+numCond))
    cdfplot(summaryStr.measureZ)
    hold on
    xValues = linspace(min(summaryStr.measureZ), max(summaryStr.measureZ));
    plot(xValues, normcdf(xValues, 0, 1), 'r-')
    legend('Empirical CDF','Standard Normal CDF','Location','best') 
end
suptitle({pA.pAnalysis, varName})

dirs.DistributionFigureFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis varName suffix 'DistributionPlot.jpg']);
export_fig(dirs.DistributionFigureFile)
end