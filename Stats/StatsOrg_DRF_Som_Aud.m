function StatsOrg_DRF_Som_Aud(dirs, StatTableSomVF, StatTableSomMN, StatTableAud)
% q4: Do participants show similar compensatory respones when only Auditory
% feedback is perturbed? 

pA.pAnalysis = 'DRF_Som_Aud';

dirs.behavioralResultTable = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'BehavioralResultTable.xlsx']);

% Question 4 %%%
% Currently Expecting Fewer 'Observations' from SomVF/SomMN compared to Aud
I = ismember(StatTableAud.SubjID, StatTableSomVF.SubjID) == 0;
StatTableAudLs = StatTableAud;
StatTableAudLs(I,:) = [];

respPer_SomVF = StatTableSomVF.RespPer;
respPer_SomMN = StatTableSomMN.RespPer;
respPer_AudLs = StatTableAudLs.RespPer;

respPerCol = [respPer_SomVF, respPer_SomMN, respPer_AudLs];
pA.condName = {'SomPert Not Masked', 'SomPert Masked', 'AudPert'};
pA.pubCondName = {'SOM Pert Not Masked',...
                  'SOM Pert Masked',...
                  'AUD Pert'};
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

measureSummaryStrs        = [];
summaryVarTableAcrossCond = table();
for i = 1:pA.numCond
    % Identify the Variable and Condition
    curCond   = pA.condName{i};
    measure   = respPerCol(:, i);
    
    measureVar.varName   = 'RespPer';
    measureVar.condition = curCond;
    measureVar.units     = 'cents';

    % Perform Standard Sumamry Stats
    summaryStat = MeasureSummaryStats(dirs, pA, measureVar, measure, lambdas(i));
    
    if ApplyTrans == 1
        usedLambda  = lambdas(2);
        summaryStat = summaryStat.performSpecificBoxCoxTrans(usedLambda);
    end
    
    % Describe the normality
    summaryStat = summaryStat.testNormality();

    % Concatenate the Summary Stat Arrays across condition
    summaryVarTableAcrossCond = [summaryVarTableAcrossCond; summaryStat.SummaryTable];

    % Concatenate the Structure for Histogram and Transformed Values
    measureSummaryStrs = cat(1, measureSummaryStrs, summaryStat.SummaryStruct);
end

drawHistograms(measureSummaryStrs, dirs, pA)

varCmp = [1 2; 1 3; 2 3];
varH   = [1.08, 1.16, 0.82];
comp = length(varCmp);

measureSummaryStrDiff = [];
for jj = 1:comp
    cond = ['DiffBetween' num2str(varCmp(jj,:))];
    % Find the difference between the two conditions and place in Struct
    measDiff = measureSummaryStrs(varCmp(jj,1)).measure - measureSummaryStrs(varCmp(jj,2)).measure;
    
    measureDiffVar.varName   = 'RespPer';
    measureDiffVar.condition = cond;
    measureDiffVar.units     = 'cents';
    summaryStatDiff = MeasureSummaryStats(dirs, pA, measureDiffVar, measDiff, 0);
    
    summaryStatDiff.SummaryStruct.vars = varCmp(jj,:);
    summaryStatDiff.SummaryStruct.h    = varH(jj);

    % Apply a Box Cox transform with the default (best) lambda
    summaryStatDiff = performSimpleBoxCoxTrans(summaryStatDiff);

    summaryStatDiff = summaryStatDiff.testNormality();       % Test Normality
    summaryStatDiff = summaryStatDiff.performTTest();        % Perform t-test
    summaryStatDiff.drawHistoBoxCombo()                      % Visualize Normality/Outliers

    measureSummaryStrDiff = cat(1, measureSummaryStrDiff, summaryStatDiff.SummaryStruct);
    
    % Save Behavioral Result Table: Values ready for inclusion in manuscript 
    writetable(summaryStatDiff.statSentTable, dirs.behavioralResultTable, 'Range', 'A7', 'WriteRowNames', 1, 'Sheet', ['RespPer' cond])
end

drawBoxPlot(dirs, pA, measureSummaryStrs, measureSummaryStrDiff)
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
axisLSize = 20;

cond      = pA.pubCondName;
numCond   = pA.numCond;

measBox = figure('Color', [1 1 1]);
plotpos = [30 0]; plotdim = [800 1000];
set(measBox, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

collData = [];
for i = 1:numCond
    collData(:,i) = measureSummaryStrs(i).measure;
end
minCol = min(min(collData)) - 5;

boxplot(collData)
ylabel('RespPer (%)')
% title({'Comparison of Response Percentages'; 'Between Experimental Conditions'})
box off

yt = get(gca, 'YTick');
axis([xlim    minCol  ceil(max(yt)*1.25)])
xt = get(gca, 'XTick');

numComp = length(summaryStrDiff);
for ii = 1:numComp
    cS = summaryStrDiff(ii);
    pVal = cS.ttestP;
    
    if pVal >= 0.001
        pVal = round(pVal, 3);
        pSentence = ['p = ' pVal];
    else
        pSentence = 'p < 0.001';
    end    
    
    barR = xt(cS.vars);
    barM = mean(barR);
    if cS.isSig == 1
        hold on
        plot(barR, [1 1]*max(yt)*(cS.h+.01), '-k', 'LineWidth', 2)
        plot(barM, max(yt)*(cS.h+.04), '*k', 'MarkerSize', 10)
    end
%     text(barM-0.30, max(yt)*(cS.h+.035), pSentence, 'FontSize', axisLSize, 'FontWeight','bold')
end

set(gca, 'XTickLabel', cond)
fix_xticklabels(gca, 0.2, {'FontSize', axisLSize, 'FontName', fontN, 'FontWeight','bold'});

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold',...
        'LineWidth', 2)

dirs.BoxPlotFigureFile = fullfile(dirs.SavResultsDir, 'Question4BoxPlot.jpg');
export_fig(dirs.BoxPlotFigureFile)
end