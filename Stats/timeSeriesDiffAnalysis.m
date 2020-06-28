function timeSeriesDiffAnalysis(dirs, pA, allSubjRes)

onsetCont      = allSubjRes.audioMf0SecCont(:,:,1);
onsetPertCond1 = allSubjRes.audioMf0SecPert{1}(:,:,1);
onsetPertCond2 = allSubjRes.audioMf0SecPert{2}(:,:,1);

time          = allSubjRes.secTime;
meanCont      = allSubjRes.audioMf0MeanCont(:,1);
meanPertCond1 = allSubjRes.audioMf0MeanPert{1}(:,1);
meanPertCond2 = allSubjRes.audioMf0MeanPert{2}(:,1);
tStep = mean(diff(time));

[numPts, numObs] = size(onsetPertCond1);

pAnalysisFix = pA.pAnalysis;
pAnalysisFix(strfind(pAnalysisFix, '_')) = '';

contTime = figure('Color', [1 1 1]);
plotpos = [30 0]; plotdim = [1200 300];
set(contTime, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

condH = []; condP = [];
for ii = 1:numPts
    [h, p] = ttest(onsetPertCond1(ii,:), onsetPertCond2(ii,:));
    condH = cat(1, condH, h);
    condP = cat(1, condP, p);
    xRange = [-tStep/2 tStep/2] + time(ii);
    if h == 1
        area(xRange, [400 400], -200, 'FaceColor', [0.8 0.8 0.8], 'EdgeAlpha', 0)
        hold on
    end
end

ax = plot(time, meanPertCond1, 'b');
hold on
bx = plot(time, meanPertCond2, 'r');
box off
axis([-0.5 1 -100 20])
title(pAnalysisFix)

legend([ax, bx], {'Perturbed notMasked Trials', 'Perturbed Masked Trials'}, 'Location', 'southwest', 'EdgeColor', 'none')

dirs.timeSeriesConds = fullfile(dirs.SavResultsDir, [pAnalysisFix 'TimeSeriesStatCompareCond.jpg']);
export_fig(dirs.timeSeriesConds)

contTime = figure('Color', [1 1 1]);
plotpos = [30 0]; plotdim = [1200 300];
set(contTime, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

condH = []; condP = [];
for ii = 1:numPts
    [h, p] = ttest(onsetPertCond1(ii,:), onsetCont(ii,:));
    condH = cat(1, condH, h);
    condP = cat(1, condP, p);
    xRange = [-tStep/2 tStep/2] + time(ii);
    if h == 1
        area(xRange, [400 400], -200, 'FaceColor', [0.8 0.8 0.8], 'EdgeAlpha', 0)
        hold on
    end
end

ax = plot(time, meanCont, 'k');
hold on
bx = plot(time, meanPertCond2, 'r');
box off
axis([-0.5 1 -100 20])
title(pAnalysisFix)

legend([ax, bx], {'Control Trials', 'Perturbed Masked Trials'}, 'Location', 'southwest', 'EdgeColor', 'none')

dirs.timeSeriesCont = fullfile(dirs.SavResultsDir, [pAnalysisFix 'TimeSeriesStatCompareCont.jpg']);
export_fig(dirs.timeSeriesCont)
end