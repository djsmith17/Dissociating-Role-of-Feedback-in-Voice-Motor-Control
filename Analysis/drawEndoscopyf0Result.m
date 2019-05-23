function drawEndoscopyf0Result()

close all
dirs = dfDirs('Dissociating-Role-of-Feedback-in-Voice-Motor-Control');

allParti = {'DRF5', 'DRF9', 'DRF12', 'DRF14', 'DRF19'};
eachTrial = [3 10 8 5 8];

for ii = 1:5
    participant = allParti{ii};
    run         = 'SFL1';
    curTrial    = eachTrial(ii);
    curRes = analyzeAndDrawResult(dirs, participant, run, curTrial);
end
end

function curRes = analyzeAndDrawResult(dirs, participant, run, curTrial)
        
dirs.ResultsParti     = fullfile(dirs.Results, participant, run);
dirs.ResultsBehavFile = fullfile(dirs.ResultsParti, [participant run 'ResultsDRF.mat']);
dirs.ResultsCodedVid  = fullfile(dirs.ResultsParti, [participant run 'EndoFrameMeasuresDJS.mat']);

load(dirs.ResultsBehavFile) % returns res
load(dirs.ResultsCodedVid) % returns CodedEndoFrameDataSet

curRes.participant      = participant;
curRes.run              = run;
curRes.curTrial         = curTrial;
curRes.trialNums        = res.allIdxFin(res.pertIdxFin);
ii = find(curRes.trialNums == curTrial);

curRes.time             = res.timef0;
curRes.sigs             = res.audioMf0TrialPert;
curRes.pertTrig         = res.pertTrigsFin;

curRes.curSig           = curRes.sigs(:,ii);
curRes.curPertTrig      = curRes.pertTrig(ii,:);
curRes.limits           = res.limitsA;

curRes.timePres         = res.timeS;
curRes.sensorP          = res.sensorPsv;
curRes.curSensorP       = curRes.sensorP(:,ii);
curRes.pressureLim      = res.limitsP;


% Video Results:
curTable = CodedEndoFrameDataSet{ii};

curRes.timeFrames = linspace(0,4,121);
pt1X = curTable.FidPt1X;
pt2X = curTable.FidPt2X;
pt1Y = curTable.FidPt1Y;
pt2Y = curTable.FidPt2Y;
distRaw = curTable.Dist;
[curRes.dist, curRes.distLim] = adjustDistMeasure(curRes.timeFrames, distRaw);

XPtsDiff = pt2X - pt1X;
YPtsDiff = pt2Y - pt1Y;

drawEndoResponses(dirs, curRes)
end

function drawEndoResponses(dirs, curRes)

plotpos = [10 0];
plotdim = [1000 500];
InterTrialf0 = figure('Color', [1 1 1]);
set(InterTrialf0, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

pertColor = [0.8 0.8 0.8];

pertAx  = [curRes.curPertTrig(1), curRes.curPertTrig(2)];
pertAy  = [600 600];

% f0 Results
subplot(2,1,1)
yyaxis right
plot(curRes.timePres, curRes.curSensorP, '--k', 'LineWidth', 1.5)
ylabel('Pressure (psi)')
axis(curRes.pressureLim);
set(gca,'FontSize', 14,...
        'FontWeight','bold')
yyaxis left

pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'EdgeColor', pertColor);
hold on    

plot(curRes.time, curRes.curSig, 'b', 'LineWidth', 2)
ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')
title({[curRes.participant ' ' curRes.run],[ ' Trial ' num2str(curRes.curTrial)]}, 'FontSize', 18, 'FontWeight', 'bold')
axis(curRes.limits); box off
lgd = legend(pA, 'Perturbation Period');
lgd.EdgeColor = 'none';

set(gca,'FontSize', 14,...
        'FontWeight','bold')

subplot(2,1,2)
pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'EdgeColor', pertColor);
hold on
plot([-1 5], [0 0], '--k')
ptDistLn = plot(curRes.timeFrames, curRes.dist, 'b');
axis([0 4 curRes.distLim]); box off
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold')
ylabel('Distance (Pixels)')
title('Euclidian Distance Between Points')

set(gca,'FontSize', 14,...
        'FontWeight','bold')

fileName = 'EndoDistanceMeasure.png';
savedFigFile = fullfile(dirs.ResultsParti, [curRes.participant curRes.run fileName]);
export_fig(savedFigFile)
end

function [dist, distLim] = adjustDistMeasure(time, distRaw)

indPrePert = time < 1.0;

distSmooth = smooth(distRaw, 3);
distBase   = mean(distSmooth(indPrePert));

dist = distSmooth - distBase;
distMax = max(dist) + 10;
distMin = min(dist) - 10;
distLim = [distMin distMax];
end