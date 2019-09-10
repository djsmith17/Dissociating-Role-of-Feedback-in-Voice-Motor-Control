function drawEndoscopyf0Result()
% drawEndoscopyf0Result() is the function that organizes the coded
% endoscopy data set. 
%
% This function calls the following subfunctions
% -analyseAndDrawResult()
% -iterateOnAnalysisSteps()
% -drawEndoResponses()

close all
dirs = dfDirs('Dissociating-Role-of-Feedback-in-Voice-Motor-Control');
dirs.PooledResultsDir = fullfile(dirs.Results, 'Pooled Analyses', 'DRF_Endo');

eAn.allParti = {'DRF5', 'DRF9', 'DRF12', 'DRF14', 'DRF19'};
eAn.numParti = length(eAn.allParti);
eAn.runs     = 'SFL1';
eAn.coder    = 'RF';

% Initalize the Pooled Endoscopy Results
eAnPool = initPooledEndoscopyResults();
for ii = 1:eAn.numParti
    [curRes, dMeasObj, dMeasObj2, dMeasObj3] = analyzeAndDrawResult(dirs, eAn.allParti{ii}, eAn.runs, eAn.coder);
    eAnPool = iterPooledEndoscopyResults(eAnPool, curRes, dMeasObj, dMeasObj2, dMeasObj3);
end
close all

eAnPool = meanPooledEndoscopyResults(eAnPool);

% Create All Subj object based on parameters from last single subj obj
distObj1AllSubj = dMeasObj;

distObj1AllSubj.curSess       = 'Mean Participant Distance Change';
distObj1AllSubj.numTrial      = eAn.numParti;
distObj1AllSubj.iterationType = 'Participants';
distObj1AllSubj.legendCurves  = [];
distObj1AllSubj.legendLabels  = {};

distObj2AllSubj = distObj1AllSubj;
distObj3AllSubj = distObj1AllSubj;

distObj1AllSubj.sigsSec = eAnPool.allSubjMeanSecs;
distObj2AllSubj.sigsSec = eAnPool.allSubjMeanSecs2;
distObj3AllSubj.sigsSec = eAnPool.allSubjMeanSecs3;

% Mean the sectioned trials
eAnPool.sigsSecf0M       = distObj1AllSubj.meanData(eAnPool.sigsSecf0);
eAnPool.codedSensorPSecM = distObj1AllSubj.meanData(eAnPool.codedSensorPSec);
distObj1AllSubj.sigsSecM = distObj1AllSubj.meanData(distObj1AllSubj.sigsSec);
distObj2AllSubj.sigsSecM = distObj2AllSubj.meanData(distObj2AllSubj.sigsSec);
distObj3AllSubj.sigsSecM = distObj3AllSubj.meanData(distObj3AllSubj.sigsSec);

eAnPool.lagTimeM  = mean(eAnPool.lagTimes, 1);
eAnPool.riseTimeM = mean(eAnPool.riseTimes, 1);

% Identify the bounds for these data
distObj1AllSubj = distObj1AllSubj.identifyBounds;
distObj1AllSubj = distObj1AllSubj.drawSigsSecM;

% Append additional lines to the plot
distObj1AllSubj = distObj1AllSubj.appendFigure(distObj2AllSubj.sigsSecM, 2);
distObj1AllSubj = distObj1AllSubj.appendFigure(distObj3AllSubj.sigsSecM, 3);

distObj1AllSubj.saveSigsSecMFig(dirs.PooledResultsDir)

% Mean Lines
distObj1AllSubj.iterationType = 'Participantsx3Lines';
distObj1AllSubj.legendCurves = [];
distObj1AllSubj.legendLabels = {};

sigsSecLinesOn = [distObj1AllSubj.sigsSecM.ON.mean,...
                  distObj2AllSubj.sigsSecM.ON.mean,...
                  distObj3AllSubj.sigsSecM.ON.mean];
              
sigsSecLinesOf = [distObj1AllSubj.sigsSecM.OF.mean,...
                  distObj2AllSubj.sigsSecM.OF.mean,...
                  distObj3AllSubj.sigsSecM.OF.mean];
              
sigsSecLines = sigsSecLinesOn;
sigsSecLines = cat(3, sigsSecLines, sigsSecLinesOf);
distObj1AllSubj.sigsSecM = distObj1AllSubj.meanData(sigsSecLines);

distObj1AllSubj = distObj1AllSubj.identifyBounds;

% All three lines collapsed: Onset Figure
stimTraceProp.timeP          = eAnPool.timePresSec;
stimTraceProp.meanPres       = eAnPool.codedSensorPSecM;
stimTraceProp.meanOnsetLag   = eAnPool.lagTimeM(1);
stimTraceProp.meanOnsetRise  = eAnPool.riseTimeM(1);
stimTraceProp.meanOffsetLag  = eAnPool.lagTimeM(2);
stimTraceProp.meanOffsetRise = eAnPool.riseTimeM(2);

distObj1AllSubj.dataColor1 = [0 78 0]/255;
distObj1AllSubj.sigsSecMLims = [-0.5 1.0 0 15.8];
distObj1AllSubj = distObj1AllSubj.drawSigsSecM_Onset(1, stimTraceProp);
distObj1AllSubj.sigsMeanFigTitle = [distObj1AllSubj.curSess '_InterTrialMeanLineOnset' distObj1AllSubj.coder '.jpg'];
distObj1AllSubj.appendFigureDynamics_Onset(eAnPool.timef0Sec, eAnPool.sigsSecf0M)

distObj1AllSubj.saveSigsSecMFig(dirs.PooledResultsDir)

distObj1AllSubj.legendCurves = [];
distObj1AllSubj.legendLabels = {};
distObj1AllSubj.sigsSecMLims = [-0.5 1.0 0 15.8];
distObj1AllSubj = distObj1AllSubj.drawSigsSecM_Offset(1, stimTraceProp);
distObj1AllSubj.sigsMeanFigTitle = [distObj1AllSubj.curSess '_InterTrialMeanLineOffset' distObj1AllSubj.coder '.jpg'];
distObj1AllSubj.appendFigureDynamics_Offset(eAnPool.timef0Sec, eAnPool.sigsSecf0M)

distObj1AllSubj.saveSigsSecMFig(dirs.PooledResultsDir)
end

function eAnPool = initPooledEndoscopyResults()
% A place to organize the extra pooled vars that need to be accounted for

eAnPool.timef0Sec       = [];
eAnPool.sigsSecf0       = [];
eAnPool.sigsSecf0M      = [];

eAnPool.timePresSec = [];
eAnPool.codedSensorPSec  = [];
eAnPool.codedSensorPSecM = [];

eAnPool.lagTimes = [];
eAnPool.lagTimeM = [];
eAnPool.riseTimes = [];
eAnPool.riseTimeM = [];

eAnPool.meanSecsOn      = [];
eAnPool.meanSecsOf      = [];
eAnPool.allSubjMeanSecs = [];

eAnPool.meanSecsOn2      = [];
eAnPool.meanSecsOf2      = [];
eAnPool.allSubjMeanSecs2 = [];

eAnPool.meanSecsOn3      = [];
eAnPool.meanSecsOf3      = [];
eAnPool.allSubjMeanSecs3 = [];
end

function eAnPool = iterPooledEndoscopyResults(eAnPool, curRes, dMeasObj, dMeasObj2, dMeasObj3)

eAnPool.timef0Sec = curRes.timef0Sec;
eAnPool.sigsSecf0 = cat(2, eAnPool.sigsSecf0, curRes.codedSigsSecM);

eAnPool.timePresSec = curRes.timePresSec;
eAnPool.codedSensorPSec = cat(2, eAnPool.codedSensorPSec, curRes.codedSensorPSecM);

eAnPool.lagTimes  = cat(1, eAnPool.lagTimes, curRes.codedLagTimesM);
eAnPool.riseTimes = cat(1, eAnPool.riseTimes, curRes.codedRiseTimesM);

eAnPool.meanSecsOn = cat(2, eAnPool.meanSecsOn, dMeasObj.sigsSecM.ON.mean);
eAnPool.meanSecsOf = cat(2, eAnPool.meanSecsOf, dMeasObj.sigsSecM.OF.mean);

if ~isempty(dMeasObj2)
    eAnPool.meanSecsOn2 = cat(2, eAnPool.meanSecsOn2, dMeasObj2.sigsSecM.ON.mean);
    eAnPool.meanSecsOf2 = cat(2, eAnPool.meanSecsOf2, dMeasObj2.sigsSecM.OF.mean);

    eAnPool.meanSecsOn3 = cat(2, eAnPool.meanSecsOn3, dMeasObj3.sigsSecM.ON.mean);
    eAnPool.meanSecsOf3 = cat(2, eAnPool.meanSecsOf3, dMeasObj3.sigsSecM.OF.mean);
end
end

function eAnPool = meanPooledEndoscopyResults(eAnPool)

eAnPool.allSubjMeanSecs = cat(3, eAnPool.allSubjMeanSecs, eAnPool.meanSecsOn);
eAnPool.allSubjMeanSecs = cat(3, eAnPool.allSubjMeanSecs, eAnPool.meanSecsOf);

eAnPool.allSubjMeanSecs2 = cat(3, eAnPool.allSubjMeanSecs2, eAnPool.meanSecsOn2);
eAnPool.allSubjMeanSecs2 = cat(3, eAnPool.allSubjMeanSecs2, eAnPool.meanSecsOf2);

eAnPool.allSubjMeanSecs3 = cat(3, eAnPool.allSubjMeanSecs3, eAnPool.meanSecsOn3);
eAnPool.allSubjMeanSecs3 = cat(3, eAnPool.allSubjMeanSecs3, eAnPool.meanSecsOf3);
end

function [curRes, dMeasObj, dMeasObj2, dMeasObj3] = analyzeAndDrawResult(dirs, participant, run, coder)

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
curRes.timef0Sec        = res.secTime;
curRes.sigsf0Sec        = res.audioMf0SecPert;
curRes.pertTrig         = res.pertTrigsFin;
curRes.limits           = res.limitsA;

% Pressure Sensor signals
curRes.timePres         = res.timeS;
curRes.sensorP          = res.sensorPsv;
curRes.pressureLim      = res.limitsP;

% Pressure Sensor dynamics
curRes.presSD = res.presSDsv;

% How many trials were coded?
curRes.numTrial         = res.numPertTrialsFin;

% Video Properties
[numFrame, ~] = size(CodedEndoFrameDataSet{1});
curRes.timeFrames = linspace(0, curRes.time(end), numFrame);

ctIdx = [];
curRes.codedDist  = []; 
curRes.codedDist2 = []; 
curRes.codedDist3 = [];
for ii = 1:curRes.numTrial
    curTable = CodedEndoFrameDataSet{ii};
    firstVal = curTable.FidPt1X;
    if firstVal ~= 0 % Was this trial coded at all? This check needs to be improved in future
        ctIdx = [ctIdx ii];
       
        curRes.codedDist     = cat(2, curRes.codedDist, curTable.Dist);
        
        if ismember('Dist2', curTable.Properties.VariableNames)
            curRes.codedDist2    = cat(2, curRes.codedDist2, curTable.Dist2);
            curRes.codedDist3    = cat(2, curRes.codedDist3, curTable.Dist3);
        end
    end    
end

curRes.codedTrialNum      = curRes.trialNums(ctIdx);
curRes.codedPertTrig      = curRes.pertTrig(ctIdx,:);
curRes.codedSigs          = curRes.sigs(:,ctIdx);
curRes.codedSigsSec       = curRes.sigsf0Sec(:, ctIdx, :);
curRes.codedSensorP       = curRes.sensorP(:, ctIdx);
curRes.timePresSec        = curRes.presSD.timeSec;
curRes.codedSensorPSec    = curRes.presSD.sensorSec(:, ctIdx, :);
curRes.codedSensorTrigTSt = curRes.presSD.TrigTime(ctIdx,:);
curRes.codedSensorTrigTSp = curRes.presSD.TrigTime(ctIdx,:) + curRes.presSD.riseTimes(ctIdx,:);

curRes.codedlagTimes  = curRes.presSD.lagTimes(ctIdx,:)- curRes.presSD.lagTimes(ctIdx,:);
curRes.codedriseTimes = curRes.presSD.riseTimes(ctIdx,:);

% Mean the behavioral stuff
curRes.codedSigsSecM    = mean(curRes.codedSigsSec, 2);
curRes.codedSensorPSecM = mean(curRes.codedSensorPSec, 2);
curRes.codedLagTimesM   = mean(curRes.codedlagTimes, 1);
curRes.codedRiseTimesM  = mean(curRes.codedriseTimes, 1);

% Set up the sectioned data object
dataInfo.curSess = curRes.curSess;
dataInfo.sigType = 'Euclidian Distance';
dataInfo.units   = 'pixels';
dataInfo.coder   = coder;
dataInfo.itrType = 'Trials';

dMeasObj = iterateOnAnalysisSteps(curRes.timeFrames, curRes.codedDist, [curRes.codedPertTrig+ curRes.presSD.lagTimes(ctIdx,:)], dataInfo);

% Draw the mean-trial Onset and Offset traces
stimWindowProp.meanOnsetLag  = curRes.presSD.lagTimeM(1)/1000;
stimWindowProp.meanOnsetRise = curRes.presSD.riseTimeM(1)/1000;
stimWindowProp.meanOffsetLag  = curRes.presSD.lagTimeM(2)/1000;
stimWindowProp.meanOffsetRise = curRes.presSD.riseTimeM(2)/1000;
dMeasObj = dMeasObj.drawSigsSecM(stimWindowProp);

if ismember('Dist2', curTable.Properties.VariableNames)
    dMeasObj2 = iterateOnAnalysisSteps(curRes.timeFrames, curRes.codedDist2, [curRes.codedPertTrig+ curRes.presSD.lagTimes(ctIdx,:)], dataInfo);
    dMeasObj  = dMeasObj.appendFigure(dMeasObj2.sigsSecM, 2);
    
    dMeasObj3 = iterateOnAnalysisSteps(curRes.timeFrames, curRes.codedDist3, [curRes.codedPertTrig+ curRes.presSD.lagTimes(ctIdx,:)], dataInfo);
    dMeasObj  = dMeasObj.appendFigure(dMeasObj3.sigsSecM, 3);
else
    dMeasObj2 = [];
    dMeasObj3 = [];
end
dMeasObj.saveSigsSecMFig(dirs.ResultsParti)

% Draw the individual trial f0 vs inflation trace
% drawEndoResponses(dirs, curRes, dMeasObj, 3)
end

function dMeasObj = iterateOnAnalysisSteps(timeFrames, codedDist, codedPertTrig, dataInfo)
% dMeasObj = iterateOnAnalysisSteps(timeFrames, codedDist, codedPertTrig, dataInfo)
% sets up the necessary class methods that need to be run for this analysis
%

fs = 30;

% Create the object that handles signal sectioning
dMeasObj = dfSectionDataOrg(timeFrames, codedDist, codedPertTrig, fs, dataInfo);

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
end

function drawEndoResponses(dirs, curRes, dMeasObj, ii)

plotpos = [10 0];
plotdim = [1200 700];
InterTrialf0 = figure('Color', [1 1 1]);
set(InterTrialf0, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

% Perturbation period
pertColor = [0.8 0.8 0.8];
pertAx = [curRes.codedPertTrig(ii, 1), curRes.codedPertTrig(ii, 2)];
pertAy = [600 600];

% Inflation period
InflColor = [0 1 0];
InflAx = [curRes.codedSensorTrigTSt(ii, 1) curRes.codedSensorTrigTSp(ii, 1)];
InflAy = [600 600];

% Deflation period
DeflColor = [1 0 0];
DeflAx = [curRes.codedSensorTrigTSt(ii, 2) curRes.codedSensorTrigTSp(ii, 2)];
DeflAy = [600 600];

% pressure trace
subplot(2,1,1)
pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
hold on 
InflA = area(InflAx, InflAy, -600, 'FaceColor', InflColor, 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
hold on
DeflA = area(DeflAx, DeflAy, -600, 'FaceColor', DeflColor, 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
hold on

% f0 trace
plot(curRes.time, curRes.codedSigs(:,ii), 'b', 'LineWidth', 2)
ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')
title({[curRes.participant ' ' curRes.run],[ ' Trial ' num2str(curRes.codedTrialNum(ii))]}, 'FontSize', 18, 'FontWeight', 'bold')
axis(curRes.limits); box off

set(gca,'FontSize', 14,...
        'FontWeight','bold')

yyaxis right
plot(curRes.timePres, curRes.codedSensorP(:, ii), '--k', 'LineWidth', 1.5)
ylabel('Pressure (psi)', 'Color', 'k')
axis(curRes.pressureLim);
set(gca,'FontSize', 14,...
        'FontWeight','bold')

% Distance Results
subplot(2,1,2)
pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
hold on 
InflA = area(InflAx, InflAy, -600, 'FaceColor', InflColor, 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
hold on
DeflA = area(DeflAx, DeflAy, -600, 'FaceColor', DeflColor, 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
hold on

plot([-1 5], [0 0], '--r')
ptDistLn = plot(dMeasObj.time, dMeasObj.sigs(:, ii), 'b', 'LineWidth', 2);
axis(dMeasObj.sigsLims); box off
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold')
ylabel('Distance (Pixels)')
title('Euclidian Distance Between Points')

lgd = legend([pA InflA DeflA], {'Perturbation Period', 'Inflation Period', 'Deflation Period'});
lgd.EdgeColor = 'none';
lgd.Location = 'SouthWest';

set(gca,'FontSize', 14,...
        'FontWeight','bold')

fileName = ['EndoDistanceMeasure' curRes.coder '.jpg'];
savedFigFile = fullfile(dirs.ResultsParti, [curRes.participant curRes.run fileName]);
export_fig(savedFigFile)
end