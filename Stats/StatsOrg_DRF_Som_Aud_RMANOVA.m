function StatsOrg_DRF_Som_Aud_RMANOVA(dirs, StatTableSomVF, StatTableSomMN, StatTableAud, statTableJND)
% q4: Do participants show similar compensatory respones when only Auditory
% feedback is perturbed? 

pA.pAnalysis = 'DRF_Som_Aud';
pA.cond      = {'SomPert Not Masked', 'SomPert Masked', 'AudPert'};
pA.pubCond   = {'Laryngeal perturbation without auditory masking',...
                'Laryngeal perturbation with auditory masking',...
                'Auditory perturbation'};
pA.numCond   = length(pA.cond);

meas    = {'StimMag', 'RespMag', 'RespPer', 'tAtMin', 'f0'};
measPub = {'Stimulus Magnitude', 'Response Magnitude', 'Response Percentage', 'Response Latency', 'Fundamental Frequency'};
mUnits  = {'cents', 'cents', '%', 's', 'Hz'};
numMeas = length(meas);

measToBeTransformedBoxCox = {'StimMag', 'RespMag'};
measLambdasToUse          = [2 1 0 0 0];

meas4NonParametric = {'tAtMin', 'f0'};

textFileName = fullfile(dirs.SavResultsDir, 'rmANOVAStats.txt');
fid = fopen(textFileName, 'w');

dirs.behavioralResultTable = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'BehavioralResultTable.csv']);

I = ismember(StatTableAud.SubjID, StatTableSomVF.SubjID) == 0;
StatTableAudLs = StatTableAud;
StatTableAudLs(I,:) = [];
statTableJNDLs = statTableJND;
statTableJNDLs(I,:) = [];
JND = statTableJNDLs.JNDScoreMean;

ApplyTrans = 1;
for k = 1:numMeas
    
    [curStatTable, cond_table] = organizeVarByCond(StatTableSomVF, StatTableSomMN, StatTableAudLs, statTableJNDLs, meas{k}, pA.cond);
    
    lambdas = [];
    if ApplyTrans == 1
        for i = 1:pA.numCond
            % Identify the Variable and Condition
            curCond   = cond_table{i};
            measure   = curStatTable.(curCond);

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
        curCond   = cond_table{i};
        measure   = curStatTable.(curCond);

        measureVar.varName    = meas{k};
        measureVar.varNamePub = measPub{k};
        measureVar.condition  = curCond;
        measureVar.units      = mUnits{k};

        % Perform Standard Sumamry Stats
        summaryStat = MeasureSummaryStats(dirs, pA, measureVar, measure, lambdas(i));

        if ismember(meas{k}, measToBeTransformedBoxCox) && ApplyTrans == 1
            usedLambda  = lambdas(measLambdasToUse(k));
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

    fprintf(fid, '### %s ###\n\n', meas{k});
    
    if ismember(meas{k}, meas4NonParametric)       
        isSig = performNonParametricTesting(fid, measureSummaryStrs);
    else
        isSig = performParametricTesting(fid, measureSummaryStrs);
    end

    varCmp = [1 2; 1 3; 2 3];
    varH   = [1.08, 1.16, 1.02;...
              1.08, 1.16, 0.82;...
              1.08, 1.16, 0.82;...
              1.08, 1.16, 1.02;...
              1.08, 1.16, 0.82];
    comp = length(varCmp);

    bigTable = table('size', [3 2], 'VariableNames', {'Comparisons', 'StatSentence'}, 'VariableTypes', {'string', 'string'});
    measureSummaryStrDiff = [];
    for jj = 1:comp
        cond = ['DiffBetween' num2str(varCmp(jj,:))];
        % Find the difference between the two conditions and place in Struct
        measDiff = measureSummaryStrs(varCmp(jj,1)).measureT - measureSummaryStrs(varCmp(jj,2)).measureT;

        measureDiffVar.varName    = meas{k};
        measureDiffVar.varNamePub = measPub{k};
        measureDiffVar.condition  = cond;
        measureDiffVar.units      = mUnits{k};
        summaryStatDiff = MeasureSummaryStats(dirs, pA, measureDiffVar, measDiff, 0);
        summaryStatDiff.SummaryStruct.suffix = sprintf('DiffBetw%s%s', pA.cond{varCmp(jj,1)}, pA.cond{varCmp(jj,2)});
        summaryStatDiff.SummaryStruct.go4PostHoc = isSig;

        summaryStatDiff.SummaryStruct.vars = varCmp(jj,:);
        summaryStatDiff.SummaryStruct.h    = varH(k, jj);

        summaryStatDiff = summaryStatDiff.testNormality();       % Test Normality
        
        if summaryStatDiff.SummaryStruct.swH == 1
            summaryStatDiff = summaryStatDiff.performSimpleBoxCoxTrans();
            summaryStatDiff = summaryStatDiff.testNormality();       % Test Normality
        end
        
        if ismember(meas{k}, meas4NonParametric)
            summaryStatDiff = summaryStatDiff.performWilcoxonRankTest(pA.cond{varCmp(jj,1)}, pA.cond{varCmp(jj,2)});        % Perform t-test
            fprintf(fid, 'The r effect size for the below Wilcoxon signed rank test is %0.3f\n', (summaryStatDiff.SummaryStruct.zstat/sqrt(summaryStatDiff.SummaryStruct.numObvs)));
        else
            summaryStatDiff = summaryStatDiff.performTTest(pA.cond{varCmp(jj,1)}, pA.cond{varCmp(jj,2)});        % Perform t-test
        end
        
        summaryStatDiff.drawHistoBoxCombo()                      % Visualize Normality/Outliers
        measureSummaryStrDiff = cat(1, measureSummaryStrDiff, summaryStatDiff.SummaryStruct);

        bigTable.Comparisons(jj)  = cond;
        bigTable.StatSentence(jj) = summaryStatDiff.statSentTable.StatSentence;
        
        fprintf(fid, summaryStatDiff.statSentTable.StatSentence{1});                                                                                                                                                            
    end

    drawBoxPlot(dirs, pA, measureSummaryStrs, measureSummaryStrDiff)
    fprintf(fid, '\n');
end

fclose(fid);
end

function [curStatTable, cond_Table] = organizeVarByCond(StatTableSomVF, StatTableSomMN, StatTableAud, StatTableJND, meas, cond)

cond_Table = matlab.lang.makeValidName(cond); % Valid idenitifer strings
condSubj = ['SubjID' cond_Table 'JND'];

subjID     = StatTableSomVF{:, 'SubjID'};
VFmeasure  = StatTableSomVF{:, meas};
MNmeasure  = StatTableSomMN{:, meas};
Audmeasure = StatTableAud{:, meas};
Acuity     = StatTableJND{:, 'JNDScoreMean'};

curStatTable = table(subjID, VFmeasure, MNmeasure, Audmeasure, Acuity, 'VariableNames', condSubj);
end

function isSig = performParametricTesting(fid, measureSummaryStrs)

measTable       = table();
measTable.SomVF = measureSummaryStrs(1).measureT;
measTable.SomMN = measureSummaryStrs(2).measureT;
measTable.Aud   = measureSummaryStrs(3).measureT;
expTableTrans   = table({'SomVF' 'SomMN', 'Aud'}','VariableNames',{'Conditions'});

% Create a model in order to test for sphericity
measFit = fitrm(measTable, 'SomVF-Aud~1', 'WithinDesign', expTableTrans, 'WithinModel', 'separatemeans');
measSph = mauchly(measFit);
Psphericty = measSph.pValue;

% From testing of normality, and applying transforms, we can say that the
% three measures meet the assumptions to perform a 
rAnovaRes = ranova(measFit);

if Psphericty >= 0.05
    Panova = rAnovaRes.pValue(1);
    adjustNote = 'We will use the un-adjusted p-Value\n';
    sphercityNote = ' not';
else
    Panova = rAnovaRes.pValueGG(1);
    adjustNote = 'We will use the Greenhouse-Geisser adjusted p-Value\n';
    sphercityNote = '';
end

if Panova < 0.05/4
    sigNote = '';
    isSig   = 1;
else
    sigNote = ' not';
    isSig   = 0;
end

meas = measureSummaryStrs(1).varName;

fprintf(fid, 'The Mauchly test indicates that the assumption of sphericity has%s been violated (X2(%d) = %0.3f p = %0.3f)\n', sphercityNote, measSph.DF, measSph.ChiStat, measSph.pValue);
fprintf(fid, adjustNote);
fprintf(fid, 'The results of the rmANOVA indicate that %s was%s significantly different between the experimental conditions (F(%d, %d) = %0.2f, p = %0.6f).\n', meas, sigNote, rAnovaRes.DF(1), rAnovaRes.DF(2), rAnovaRes.F(1), Panova);

%%% Effect Size %%%
partialEta = rAnovaRes.SumSq(1)/(rAnovaRes.SumSq(1) + rAnovaRes.SumSq(2));
fprintf(fid, 'The partial eta for this rmANOVA was calculated to be %0.3f\n\n', partialEta);
end

function isSig = performNonParametricTesting(fid, measureSummaryStrs)

meas = measureSummaryStrs(1).varName;

measTable       = table();
measTable.SomVF = measureSummaryStrs(1).measureT;
measTable.SomMN = measureSummaryStrs(2).measureT;
measTable.Aud   = measureSummaryStrs(3).measureT;

[P, Ftable, ~] = friedman(measTable{:,:}, 1, 'off');

if P < 0.05/3
    sigNote = '';
    isSig   = 1;
else
    sigNote = ' not';
    isSig   = 0;
end

statSentence = fprintf(fid, 'Statistical analyses revealed that there was%s a significant effect of experimental condition on %s (chiSq(%d) = %0.2f, p = %0.6f)\n',...
                       sigNote,...
                       meas,...
                       Ftable{2,3},...
                       Ftable{2,2},...
                       P);

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
    
    step = (summaryStr.max-summaryStr.min)*0.1;
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
    
    histoTitle = [cond{ii} ' (H = ' num2str(summaryStr.swH) ', WStat = ' num2str(summaryStr.swTest) ')'];
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

cond      = pA.pubCond;
numCond   = pA.numCond;

measBox = figure('Color', [1 1 1]);
plotpos = [30 50]; plotdim = [650 750];
set(measBox, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

collData = [];
for i = 1:numCond
    collData(:,i) = measureSummaryStrs(i).measure;
end

boxplot(collData)
ylabel([measureSummaryStrs(1).varNamePub ' (' measureSummaryStrs(1).units ')'])
% title({'Comparison of Response Percentages'; 'Between Experimental Conditions'})
box off

yt = get(gca, 'YTick');
minYLim = min(min(collData))*0.75;
maxYLim = max(max(collData))*1.25;

axis([xlim minYLim  350])
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
    if cS.isSig == 1 && cS.go4PostHoc == 1
        hold on
        plot(barR, [1 1]*max(yt)*(cS.h+.01), '-k', 'LineWidth', 2)
        plot(barM, max(yt)*(cS.h+.04), '*k', 'MarkerSize', 10)
    end
%     text(barM-0.30, max(yt)*(cS.h+.035), pSentence, 'FontSize', axisLSize, 'FontWeight','bold')
end

set(gca, 'XTickLabel', cond)
fix_xticklabels(gca, 0.17, {'FontSize', 17, 'FontName', fontN, 'FontWeight','bold'});

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold',...
        'LineWidth', 2)

dirs.BoxPlotFigureFile = fullfile(dirs.SavResultsDir, [measureSummaryStrs(1).varName 'TTestBoxPlot.jpg']);
export_fig(dirs.BoxPlotFigureFile, '-r300')
end