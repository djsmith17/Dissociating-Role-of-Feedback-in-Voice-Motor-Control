function drawEndoscopyf0Result()

close all
dirs = dfDirs('Dissociating-Role-of-Feedback-in-Voice-Motor-Control');
dirs.PooledResultsDir = fullfile(dirs.Results, 'Pooled Analyses', 'DRF_Endo');

allParti = {'DRF5', 'DRF9', 'DRF12', 'DRF14', 'DRF19'};
eachTrial = [3 10 8 5 8];

meanSecsOn = [];
meanSecsOf = [];
allSubjMeanSecs = [];
for ii = 1:5
    participant = allParti{ii};
    run         = 'SFL1';
    curTrial    = eachTrial(ii);
    
    [~, dMeasObj] = analyzeAndDrawResult(dirs, participant, run, curTrial);

    meanSecsOn = cat(2, meanSecsOn, dMeasObj.sigsSecM(:,1));
    meanSecsOf = cat(2, meanSecsOf, dMeasObj.sigsSecM(:,3));
end
close all

allSubjMeanSecs = cat(3, allSubjMeanSecs, meanSecsOn);
allSubjMeanSecs = cat(3, allSubjMeanSecs, meanSecsOf);
[~, numSubj, ~] = size(allSubjMeanSecs);

% Create All Subj object based on parameters from last single subj obj
distObjAllSubj = dMeasObj;

distObjAllSubj.curSess = 'Mean Participant Distance Change';
distObjAllSubj.sigsSec = allSubjMeanSecs;
distObjAllSubj.numTrial = numSubj;

% Perform the mean on the sectioned trials
distObjAllSubj.sigsSecM = distObjAllSubj.meanData(distObjAllSubj.sigsSec);
% Identify the bounds for these data
distObjAllSubj = distObjAllSubj.identifyBounds;

distObjAllSubj = distObjAllSubj.drawSigsSecM;
distObjAllSubj.saveSigsSecMFig(dirs.PooledResultsDir)
end

function [curRes, dMeasObj] = analyzeAndDrawResult(dirs, participant, run, curTrial)

dirs.ResultsParti     = fullfile(dirs.Results, participant, run);
dirs.ResultsBehavFile = fullfile(dirs.ResultsParti, [participant run 'ResultsDRF.mat']);
dirs.ResultsCodedVid  = fullfile(dirs.ResultsParti, [participant run 'EndoFrameMeasuresDJS.mat']);

load(dirs.ResultsBehavFile) % returns res
load(dirs.ResultsCodedVid) % returns CodedEndoFrameDataSet

curRes.participant      = participant;
curRes.run              = run;
curRes.curTrial         = curTrial;
curRes.curSess          = res.curSess;
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

timeFrames = linspace(0,4,121);
pt1X = curTable.FidPt1X;
pt2X = curTable.FidPt2X;
pt1Y = curTable.FidPt1Y;
pt2Y = curTable.FidPt2Y;
distRaw = curTable.Dist;

dataInfo.curSess = curRes.curSess;
dataInfo.sigType = 'Euclidian Distance';
dataInfo.units   = 'pixels';

% Create the object that handles signal sectioning
dMeasObj = dfSectionDataOrg(timeFrames, distRaw, curRes.curPertTrig, 30, dataInfo);

% Section the raw data for ease of baseline detection
sigsSec4Baseline  = dMeasObj.sectionData(dMeasObj.sigs);

% Identify Baseline Values
dMeasObj = dMeasObj.identifyBaselineValues(sigsSec4Baseline);

% Apply smoothing
dMeasObj.sigs = dMeasObj.applySmoothing(dMeasObj.sigs, 3);

% Apply a zero-mean to the full signals, replace the original full signal
dMeasObj.sigs = dMeasObj.applyZeroMean(dMeasObj.sigs, dMeasObj.sigsBase);

% Section the processed data
dMeasObj.sigsSec  = dMeasObj.sectionData(dMeasObj.sigs);

% Perform the mean on the sectioned trials
dMeasObj.sigsSecM = dMeasObj.meanData(dMeasObj.sigsSec);
% Identify the bounds for these data
dMeasObj = dMeasObj.identifyBounds;
% Draw the mean-trial Onset and Offset traces
dMeasObj = dMeasObj.drawSigsSecM;
dMeasObj.saveSigsSecMFig(dirs.ResultsParti)

drawEndoResponses(dirs, curRes, dMeasObj)
end

function drawEndoResponses(dirs, curRes, dMeasObj)

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

% Distance Results
subplot(2,1,2)
pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'EdgeColor', pertColor);
hold on
plot([-1 5], [0 0], '--k')
ptDistLn = plot(dMeasObj.time, dMeasObj.sigs, 'b');
axis(dMeasObj.sigsLims); box off
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold')
ylabel('Distance (Pixels)')
title('Euclidian Distance Between Points')

set(gca,'FontSize', 14,...
        'FontWeight','bold')

fileName = 'EndoDistanceMeasure.png';
savedFigFile = fullfile(dirs.ResultsParti, [curRes.participant curRes.run fileName]);
export_fig(savedFigFile)
end