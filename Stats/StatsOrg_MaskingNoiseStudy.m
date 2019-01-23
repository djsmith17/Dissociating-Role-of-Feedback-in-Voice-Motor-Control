function StatsOrg_MaskingNoiseStudy(dirs, pA, allSubjRes)

dirs.behavioralResultTable = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'BehavioralResultTable.xlsx']);

allSubjStatTable = allSubjRes.statTable;
meas = {'StimMag', 'RespMag', 'RespPer'};
units = {'cents', 'cents', '%'};
colors = ['b', 'r', 'g'];
mu     = '\mu';
sigma  = '\sigma';
nMeas = length(meas);

cond    = pA.cond;
numCond = pA.numCond;

curTestingMeas = 1;

allMeasureStats = [];
for k = curTestingMeas
    curStatTable = allSubjStatTable(:, {'AudFB', meas{k}});
    
    measDist = figure('Color', [1 1 1]);
    plotpos = [10 10]; plotdim = [1300 800];
    set(measDist, 'Position',[plotpos plotdim],'PaperPositionMode','auto')
    ha = tight_subplot(2, numCond,[0.1 0.05],[0.1 0.1],[0.05 0.05]);

    measStats = [];
    for i = 1:numCond
        curFB = strcmp(curStatTable.AudFB, cond(i));
        
        measure   = curStatTable{curFB, 2};
        
        [summaryStr, pooledVariable] = RawSummaryStats(measure);
        
        if k == 4
            summaryStr.measureT = boxcox(summaryStr.measure);
            suffix = 'Trans';
        else
            summaryStr.measureT = summaryStr.measure;
            suffix = '';
        end
                
        summaryStr.measureSkew     = round(skewness(summaryStr.measureT), 4);
        summaryStr.measureKurtosis = round(kurtosis(summaryStr.measureT), 2);
        
        summaryStr.measureZ        = zscore(summaryStr.measureT);
        [swH, swPValue, swTest]    = swtest(summaryStr.measureZ);
        
        summaryStr.swH      = swH;
        summaryStr.swPValue = round(swPValue, 3);
        summaryStr.swTest   = round(swTest, 3);
        
        pooledVariable = concateAdditionalStat(summaryStr, pooledVariable);
        
        measStat = [measureM;...
                    measureMin;...
                    measureMed;...
                    measureMax;...
                    measureSD;...
                    measureSE;...
                    measureSkew;...
                    measureKurotsis;...
                    swH;...
                    swPValue;...
                    swTest];
        measStats = cat(2, measStats, summaryStr);
        
        step = 25;
        minBound = floor(measureMin/step)*step;
        maxBound = ceil(measureMax/step)*step;
        distBin = minBound:step:maxBound;
        nBins   = length(distBin)-1;
        
        axes(ha(i))
        if k == 4
            histogram(measureTrans, 8, 'FaceColor', colors(i), 'EdgeColor', colors(i))
        else
            histogram(measureTrans, nBins, 'FaceColor', colors(i), 'EdgeColor', colors(i), 'BinLimits', [minBound maxBound])
        end
        title(cond{i})
        box off;
        
        axes(ha(i+numCond))
        cdfplot(measureZScore)
        hold on
        xValues = linspace(min(measureZScore), max(measureZScore));
        plot(xValues, normcdf(xValues, 0, 1), 'r-')
        legend('Empirical CDF','Standard Normal CDF','Location','best')        
    end
    suptitle({pA.pAnalysis, meas{k}})
    
%     measFit = fitrm(curStatTable, 'RespMag~AudFB');
%     measSph = mauchly(measFit);
    
    
    allMeasureStats = cat(2, allMeasureStats, measStats);
    dirs.DistributionFigureFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis meas{k} suffix 'DistributionPlot.jpg']);
    export_fig(dirs.DistributionFigureFile)
    
    xlswrite(dirs.behavioralResultTable, measStats, meas{k})
end

fullResultFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'AllVariableTable.xlsx']);
writetable(allSubjStatTable, fullResultFile, 'WriteVariableNames',true)
end

function [summaryStr, pooledVariable] = RawSummaryStats(measure)

numObs    = length(measure);

% Calculate the Descriptive Stats
summaryStr.measure  = measure;        % Raw Data Values
summaryStr.measureT = [];             % Transformed Data Values
summaryStr.measureZ = [];             % Z-Scored Data Values 
summaryStr.mean     = round(mean(measure), 2);
summaryStr.median   = round(median(measure), 2);
summaryStr.min      = round(min(measure), 2);
summaryStr.max      = round(max(measure), 2);
summaryStr.SD       = round(std(measure), 2);
summaryStr.SE       = round(summaryStr.SD/sqrt(numObs), 2);

pooledVariable = {summaryStr.mean;...
                  summaryStr.median;...
                  summaryStr.min;...
                  summaryStr.max;...
                  summaryStr.SD;...
                  summaryStr.SE};
                  

end

function pooledVariable = concateAdditionalStat(summaryStr, pooledVariable)

pooledVariable = cat(1, pooledVariable, summaryStr.measureSkew);
pooledVariable = cat(1, pooledVariable, summaryStr.measureKurtosis);
pooledVariable = cat(1, pooledVariable, summaryStr.swH);
pooledVariable = cat(1, pooledVariable, summaryStr.swPValue);
pooledVariable = cat(1, pooledVariable, summaryStr.swTest);
end