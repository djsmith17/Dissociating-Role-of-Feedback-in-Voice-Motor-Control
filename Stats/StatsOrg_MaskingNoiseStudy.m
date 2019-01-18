function StatsOrg_MaskingNoiseStudy(dirs, pA, allSubjRes)

dirs.behavioralResultTable = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'BehavioralResultTable.xlsx']);

allSubjStatTable = allSubjRes.statTable;
n = height(allSubjStatTable);
meas = {'StimMag', 'RespMag', 'RespPer'};
units = {'cents', 'cents', '%'};
colors = ['b', 'r', 'g'];
mu     = '\mu';
sigma  = '\sigma';
nMeas = length(meas);

cond    = pA.cond;
numCond = pA.numCond;

allMeasureStats = [];
for k = 1:nMeas
    curStatTable = allSubjStatTable(:, {'AudFB', meas{k}});
    
%     measFit = fitrm(curStatTable, 'StimMag~AudFB');
%     measSph = mauchly(measFit);
    
    measDist = figure('Color', [1 1 1]);
    plotpos = [10 10]; plotdim = [1300 800];
    set(measDist, 'Position',[plotpos plotdim],'PaperPositionMode','auto')
    ha = tight_subplot(2, numCond,[0.1 0.05],[0.1 0.1],[0.05 0.05]);

    measStats = [];
    for i = 1:numCond
        curFB = strcmp(curStatTable.AudFB, cond(i));
        
        measure   = curStatTable{curFB, 2};        
        numObs    = length(measure);
        
        if k == 1
            measureTrans = log10(measure);
        else
            measureTrans = measure;
        end
        
        measureM  = round(mean(measureTrans), 2);
        measureMed = round(median(measureTrans), 2);
        measureMin = round(min(measureTrans), 2);
        measureMax = round(max(measureTrans), 2);
        measureSD = round(std(measureTrans), 2);
        measureSE = round(measureSD/sqrt(numObs), 2);
        
        measureSkew     = round(skewness(measure), 4);
        measureKurotsis = round(kurtosis(measure), 2);
        
        measureZScore   = (measureTrans-measureM)./measureSD;
        [swH, swPValue, swTest] = swtest(measureZScore);
        
        swPValue = round(swPValue, 3);
        swTest   = round(swTest, 3);
        
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
        measStats = cat(2, measStats, measStat);
        
        step = 25;
        minBound = floor(measureMin/step)*step;
        maxBound = ceil(measureMax/step)*step;
%         distBin = minBound:step:maxBound;
%         nBins = length(distBin)-1;
        
        axes(ha(i))
        histogram(measureTrans, 8, 'FaceColor', colors(i), 'EdgeColor', colors(i))
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
    
    allMeasureStats = cat(2, allMeasureStats, measStats);
    dirs.DistributionFigureFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis meas{k} 'DistributionPlotTrans.jpg']);
    export_fig(dirs.DistributionFigureFile)
    
    xlswrite(dirs.behavioralResultTable, measStats, meas{k})
end

fullResultFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'AllVariableTable.xlsx']);
writetable(allSubjStatTable, fullResultFile, 'WriteVariableNames',true)
end