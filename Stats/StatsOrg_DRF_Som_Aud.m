function StatsOrg_DRF_Som_Aud(dirs, StatTableSomVF, StatTableSomMN, StatTableAud)
% q4: Do participants show similar compensatory respones when only Auditory
% feedback is perturbed? 

pA.pAnalysis = 'Butts';

% Question 4 %%%
% Currently Expecting Fewer 'Observations' from SomVF and SomMN
I = ismember(StatTableAud.SubjID, StatTableSomVF.SubjID) == 0;
StatTableAudLs = StatTableAud;
StatTableAudLs(I,:) = [];

respPer_SomVF = StatTableSomVF.RespPer;
respPer_SomMN = StatTableSomMN.RespPer;
respPer_AudLs = StatTableAudLs.RespPer;

respPerCol = [respPer_SomVF, respPer_SomMN, respPer_AudLs];
pA.condName = {'SomPert Not Masked', 'SomPert Masked', 'AudPert'};
[~, pA.numCond] = size(respPerCol);

ApplyTrans = 1;
lambdas = [];
if ApplyTrans == 1
    for i = 1:pA.numCond
        % Identify the Variable and Condition
        measure   = respPerCol(:,i);

        [~, lambda] = boxcox(measure + 1 - min(measure));
        lambdas = cat(1, lambdas, lambda);
    end
else
    lambdas = [0 0 0];
end

measureSummaryStrs     = [];
variableStatAcrossCond = [];
for i = 1:pA.numCond
    % Identify the Variable and Condition
    measure   = respPerCol(:, i);

    % Perform Standard Sumamry Stats
    [summaryStr, variableStat] = RawSummaryStats('RespPer', measure, lambdas(i));
    
    if ApplyTrans == 1
        usedLambda = lambdas(2);
        summaryStr.measureT   = boxcox(usedLambda, summaryStr.measure);
        summaryStr.isTrans    = 1;
        summaryStr.suffix     = 'TransVF';
        summaryStr.usedLambda = num2str(usedLambda);
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

drawHistograms(measureSummaryStrs, dirs, pA)
testNonParametric(measureSummaryStrs)

varCmp = [1 2; 1 3; 2 3];
varH   = [1.04, 1.12, 0.82];
comp = length(varCmp);

measureSummaryStrDiff = [];
for jj = 1:comp
    
    % Find the difference between the two conditions and place in Struct
    measDiff = measureSummaryStrs(varCmp(jj,1)).measure - measureSummaryStrs(varCmp(jj,2)).measure;
    [summaryStrDiff, ~] = RawSummaryStats(['Diff'], measDiff, 0);
    
    summaryStrDiff.vars = varCmp(jj,:);
    summaryStrDiff.h    = varH(jj);

    [summaryStrDiff.measureT, l] = boxcox(summaryStrDiff.measure + 1 - min(summaryStrDiff.measure));
    summaryStrDiff.isTrans    = 1;
    summaryStrDiff.suffix     = 'Trans';
    summaryStrDiff.usedLambda = num2str(round(l,2));

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

    drawHistoBoxCombo(summaryStrDiff, dirs, pA) % Visualize the Normality/Outliers of 
    
    measureSummaryStrDiff = cat(1, measureSummaryStrDiff, summaryStrDiff);
end

drawBoxPlot(dirs, pA, measureSummaryStrs, measureSummaryStrDiff)

% dirs.behavioralResultTable = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'BehavioralResultTable' summaryStr.suffix '.xlsx']);
% xlswrite(dirs.behavioralResultTable, variableStatAcrossCond, meas{k})

% fullResultFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'AllVariableTable.xlsx']);
% writetable(allSubjStatTable, fullResultFile, 'WriteVariableNames',true)
end

function [summaryStr, variableStat] = RawSummaryStats(variableName, measure, idealLambda)

numObs    = length(measure);

% Calculate the Descriptive Stats
summaryStr.varName  = variableName;
summaryStr.measure  = measure;        % Raw Data Values

summaryStr.mean     = round(mean(measure), 2);
summaryStr.median   = round(median(measure), 2);
summaryStr.min      = round(min(measure), 2);
summaryStr.max      = round(max(measure), 2);
summaryStr.SD       = round(std(measure), 2);
summaryStr.SE       = round(summaryStr.SD/sqrt(numObs), 2);

summaryStr.isTrans     = 0;       % Default is not transformed
summaryStr.measureT    = measure; % Transformed Data Values (Default is the same)
summaryStr.measureZ    = [];      % Z-Scored Data Values
summaryStr.idealLambda = idealLambda;
summaryStr.usedLambda  = 'N/A';   % Default is not transformed
summaryStr.suffix      = '';      % Default is not transformed

variableStat = {summaryStr.mean;...
                summaryStr.min;...
                summaryStr.median;...
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

function drawHistograms(measureSummaryStrs, dirs, pA)

units  = {'cents', 'cents', '%'};
colors = ['b', 'r', 'g'];
sigma  = '\sigma'; mu = '\mu';
lambda = '\lambda';

pAnalysis = pA.pAnalysis;
cond      = pA.condName;
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

function drawBoxPlot(dirs, pA, measureSummaryStrs, summaryStrDiff)

fontN = 'Arial';
axisLSize = 15;

cond      = pA.condName;
numCond   = pA.numCond;

measBox = figure('Color', [1 1 1]);
plotpos = [30 0]; plotdim = [800 1000];
set(measBox, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

collData = [];
for i = 1:numCond
    collData(:,i) = measureSummaryStrs(i).measure;
end

boxplot(collData, 'Labels', cond)
ylabel('RespPer (%)')
title('Comparison of Response Percentage between Experimental Conditions')
box off

yt = get(gca, 'YTick');
axis([xlim    0  ceil(max(yt)*1.2)])
xt = get(gca, 'XTick');

numComp = length(summaryStrDiff);
for ii = 1:numComp
    cS = summaryStrDiff(ii);
    barR = xt(cS.vars);
    barM = mean(barR);
    if cS.isSig == 1
        hold on
        plot(barR, [1 1]*max(yt)*(cS.h+.01), '-k')
        plot(barM, max(yt)*(cS.h+.02), '*k', 'MarkerSize', 10)
    end
    text(barM-0.35, max(yt)*(cS.h+.05), ['p = ' cS.fPround], 'FontSize',18)
end
set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold')

dirs.BoxPlotFigureFile = fullfile(dirs.SavResultsDir, 'Question4BoxPlot.jpg');
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

function [rAnovaRes, measSph] = testParametric(curStatTable, cond_table)

condTable = table(cond_table');

measFit = fitrm(curStatTable, 'VoiceFeedback-AC_BCMaskingNoise~1', 'WithinDesign', condTable, 'WithinModel', 'separatemeans');
measSph = mauchly(measFit);
    
rAnovaRes = ranova(measFit);
end

function [tFried] = testNonParametric(measureSummaryStrs)

measuresColl = [];
for i = 1:3
    measuresColl = cat(2, measuresColl, measureSummaryStrs(i).measureT);
end

[pFried,tFried, stats] = friedman(measuresColl);
end