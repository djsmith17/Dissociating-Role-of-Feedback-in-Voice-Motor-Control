function drawEndoscopyf0Result()

close all
dirs = dfDirs('Dissociating-Role-of-Feedback-in-Voice-Motor-Control');
dirs.PooledResultsDir = fullfile(dirs.Results, 'Pooled Analyses', 'DRF_Endo');

allParti = {'DRF5', 'DRF9', 'DRF12', 'DRF14', 'DRF19'};
eachTrial = [3 10 8 5 8];

allRes = [];
allSections = [];
for ii = [1 4 5]
    participant = allParti{ii};
    run         = 'SFL1';
    curTrial    = eachTrial(ii);
    curRes = analyzeAndDrawResult(dirs, participant, run, curTrial);
    
    if ii == 5 %This should be fixed in the future
        curRes.curPertTrig(2) = curRes.curPertTrig(2) - 1/30;
    end
    
    [curRes.secTime, curRes.secSigs] = sectionData(curRes.timeFrames, curRes.dist, curRes.curPertTrig);
    
    allSections = cat(2, allSections, curRes.secSigs);
    
    allRes = cat(1, allRes, curRes);
end
close all

poolRes.curSess = 'Mean Participant Distance Change';
poolRes.f0b     = 0;
poolRes.AudFB   = 'Auditory Masked';
poolRes.numContTrialsFin = 0;
poolRes.numPertTrialsFin = length(allRes);

poolRes.secTime = curRes.secTime;
poolRes.audioMf0MeanPert = meanAudioData(allSections);
poolRes.audioMf0MeanCont = zeros(size(poolRes.audioMf0MeanPert));

poolRes.pltName     = 'meanParticipantDistChange';
poolRes.limitsAmean = [-0.5 1.0 -10 40];
poolRes.audioDynamics = [];

poolRes.secTimeP = 0;
poolRes.sensorPAdjust = 0;
poolRes.InflDeflT = 0;

drawMeanTrialMicf0(poolRes, dirs.PooledResultsDir, 0)
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

function [secTime, secSigs] = sectionData(time, sigs, trigs)
% [secTime, secSigs] = sectionData(time, sigs, trigs) sections
% time series data around important points in time.
% 
% time:  Vector of time points (numSamp)
% sigs:  Matrix of signals to be sectioned (numSamp x numTrial)
% trigs: Onset and Offset time tiggers (numTrial x 2)
%
% secTime: Vector of time points corresponding to the sectioned window (numSampSec)
% secSigs: 3D mat of sectioned sigs (numSampSec x numTrial x event)
%          The 1st 3D layer are Onset Sections
%          The 2nd 3D later are Offset Sections

[~, numTrial] = size(sigs);
preEve  = 0.5; posEve = 1.0;

secSigs    = [];
OnsetSecs  = [];
OffsetSecs = [];
if numTrial > 0
    for ii = 1:numTrial
        OnsetT   = trigs(ii, 1); % Onset time
        OffsetT  = trigs(ii, 2); % Offset time

        OnsetTSt = round(OnsetT - preEve, 3);   % PreOnset time, rounded to nearest ms
        OnsetTSp = round(OnsetT + posEve, 3);   % PostOnset time, rounded to nearest ms
        OnsetSpan = time >= OnsetTSt & time <= OnsetTSp; % Indices corresponding to Onset period

        OffsetTSt = round(OffsetT - preEve, 3); % PreOffset time, rounded to nearest ms
        OffsetTSp = round(OffsetT + posEve, 3); % PostOffset time, rounded to nearest ms
        OffsetSpan = time >= OffsetTSt & time <= OffsetTSp; % Indices corresponding to Offset period

        OnsetSec  = sigs(OnsetSpan, ii);  % Data sectioned around Onset
        OffsetSec = sigs(OffsetSpan, ii); % Data sectioned around Offset

        OnsetSecs  = cat(2, OnsetSecs, OnsetSec);   % Sectioned signal onsets concatenated
        OffsetSecs = cat(2, OffsetSecs, OffsetSec); % Sectioned signal offsets concatenated
    end
    [numSampSec, ~] = size(OnsetSecs); % number of samples in sectioned signals
else
    numSampSec = 301;
end

secTime = linspace(-preEve, posEve, numSampSec); % time vector correspnding to the sectioned signals
secSigs(:,:,1) = OnsetSecs;  % 1st 3D layer
secSigs(:,:,2) = OffsetSecs; % 2nd 3D layer
end

function meanAudio = meanAudioData(secAudio)
% Some simple statistics on the sectioned audio data. 
% meanAudio is a vector containing the following information
% meanAudio(1) = mean Onset pitch contour
% meanAudio(2) = 95% CI of the mean Onset Pitch Contour
% meanAudio(3) = mean Offset pitch contour
% meanAudio(4) = 95% CI of the mean Offset Pitch Contour

OnsetSecs  = secAudio(:,:,1);
OffsetSecs = secAudio(:,:,2);
[~, numTrial] = size(OnsetSecs);

meanOnset  = nanmean(OnsetSecs, 2);  % across columns
meanOffset = nanmean(OffsetSecs, 2); % across columns

stdOnset   = nanstd(OnsetSecs, 0, 2);  % across columns
stdOffset  = nanstd(OffsetSecs, 0, 2); % across columns

SEMOnset   = stdOnset/sqrt(numTrial);  % Standard Error
SEMOffset  = stdOffset/sqrt(numTrial); % Standard Error

% NCIOnset   = 1.96*SEMOnset;  % 95% Confidence Interval
% NCIOffset  = 1.96*SEMOffset; % 95% Confidence Interval

meanAudio = [meanOnset SEMOnset meanOffset SEMOffset];
end