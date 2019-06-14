function drawEndoscopyf0Result()

close all
dirs = dfDirs('Dissociating-Role-of-Feedback-in-Voice-Motor-Control');
dirs.PooledResultsDir = fullfile(dirs.Results, 'Pooled Analyses', 'DRF_Endo');

allParti = {'DRF5', 'DRF9', 'DRF12', 'DRF14', 'DRF19'};
numParti = length(allParti);
coder    = 'DJS';
% eachTrial = [3 10 8 5 8];

meanSecsOn = [];
meanSecsOf = [];
allSubjMeanSecs = [];
for ii = 1:numParti
    participant = allParti{ii};
    run         = 'SFL1';
    
    [~, dMeasObj] = analyzeAndDrawResult(dirs, participant, run, coder);
    
%     drawEndoResponses(dirs, curRes, dMeasObj, eachTrial(ii))

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

function [curRes, dMeasObj] = analyzeAndDrawResult(dirs, participant, run, coder)

dirs.ResultsParti     = fullfile(dirs.Results, participant, run);
dirs.ResultsBehavFile = fullfile(dirs.ResultsParti, [participant run 'ResultsDRF.mat']);
dirs.ResultsCodedVid  = fullfile(dirs.ResultsParti, [participant run 'EndoFrameMeasures' coder '.mat']);

load(dirs.ResultsBehavFile) % returns res
load(dirs.ResultsCodedVid)  % returns CodedEndoFrameDataSet

% Participant and Run information
curRes.coder            = coder;
curRes.participant      = participant;
curRes.run              = run;
curRes.curSess          = res.curSess;
curRes.trialNums        = res.allIdxFin(res.pertIdxFin);

% General experimental parameters
curRes.time             = res.timef0;
curRes.sigs             = res.audioMf0TrialPert;
curRes.pertTrig         = res.pertTrigsFin;
curRes.limits           = res.limitsA;

curRes.timePres         = res.timeS;
curRes.sensorP          = res.sensorPsv;
curRes.pressureLim      = res.limitsP;

% How many trials were coded?
curRes.numTrial         = res.numPertTrialsFin;

% Video Properties
[numFrame, ~] = size(CodedEndoFrameDataSet{1});
curRes.timeFrames = linspace(0, curRes.time(end), numFrame);

curRes.codedPertTrig = [];
curRes.codedSigs     = [];
curRes.codedSensorP  = [];
curRes.codedDist     = [];
for ii = 1:curRes.numTrial
    curTable = CodedEndoFrameDataSet{ii};
    firstVal = curTable.FidPt1X;
    if firstVal ~= 0 % Was this trial coded at all? This check needs to be improved in future
        curRes.codedPertTrig = cat(1, curRes.codedPertTrig, curRes.pertTrig(ii,:));
        curRes.codedSigs     = cat(2, curRes.codedSigs, curRes.sigs(:,ii));
        curRes.codedSensorP  = cat(2, curRes.codedSensorP, curRes.sensorP(:,ii));
        curRes.codedDist     = cat(2, curRes.codedDist, curTable.Dist);
    end    
end

% Set up the sectioned data object
dataInfo.curSess = curRes.curSess;
dataInfo.sigType = 'Euclidian Distance';
dataInfo.units   = 'pixels';
dataInfo.coder   = coder;

% Create the object that handles signal sectioning
dMeasObj = dfSectionDataOrg(curRes.timeFrames, curRes.codedDist, curRes.codedPertTrig, 30, dataInfo);

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
end

function drawEndoResponses(dirs, curRes, dMeasObj, curTrial)

ii = find(curRes.trialNums == curTrial);

plotpos = [10 0];
plotdim = [1000 500];
InterTrialf0 = figure('Color', [1 1 1]);
set(InterTrialf0, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

pertColor = [0.8 0.8 0.8];

pertAx  = [curRes.codedPertTrig(ii, 1), curRes.codedPertTrig(ii, 2)];
pertAy  = [600 600];

% f0 Results
subplot(2,1,1)
yyaxis right
plot(curRes.timePres, curRes.codedSensorP(:, ii), '--k', 'LineWidth', 1.5)
ylabel('Pressure (psi)')
axis(curRes.pressureLim);
set(gca,'FontSize', 14,...
        'FontWeight','bold')
yyaxis left

pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'EdgeColor', pertColor);
hold on    

plot(curRes.time, curRes.curSig(:,ii), 'b', 'LineWidth', 2)
ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')
title({[curRes.participant ' ' curRes.run],[ ' Trial ' num2str(curTrial)]}, 'FontSize', 18, 'FontWeight', 'bold')
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
ptDistLn = plot(dMeasObj.time, dMeasObj.sigs(:, ii), 'b');
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