function dfIdentifyIntraCoderReliabilityScore()

close all
dirs = dfDirs('Dissociating-Role-of-Feedback-in-Voice-Motor-Control');
dirs.PooledResultsDir = fullfile(dirs.Results, 'Pooled Analyses', 'DRF_Endo');

eAn.allParti = {'DRF5', 'DRF14', 'DRF19'};
eAn.numParti = length(eAn.allParti);
eAn.runs     = 'SFL1';
eAn.coders   = {'RF', 'RF2'};

varNames = {'Participant', 'Run', 'Trial', 'Frame', 'Line', 'EuclidianDist1', 'EuclidianDist2', 'Diff'};
varTypes = {'string', 'string', 'string', 'double' 'string', 'double', 'double', 'double'};
numVar = length(varNames);

tableRow = 1;
for ii = 1:eAn.numParti
    [curRes, dMeasObj, dMeasObj2, dMeasObj3] = analyzeAndDrawResult(dirs, eAn.allParti{ii}, eAn.runs, eAn.coders)
end


end

function [curRes, dMeasObj, dMeasObj2, dMeasObj3] = analyzeAndDrawResult(dirs, participant, run, coder)

dirs.ResultsParti     = fullfile(dirs.Results, participant, run);
dirs.ResultsCodedVid1 = fullfile(dirs.ResultsParti, [participant run 'EndoFrameMeasures' coder{1} '.mat']);
dirs.ResultsCodedVid2 = fullfile(dirs.ResultsParti, [participant run 'EndoFrameMeasures' coder{2} '.mat']);

CodedSet1 = load(dirs.ResultsCodedVid1); % returns CodedEndoFrameDataSet

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

end