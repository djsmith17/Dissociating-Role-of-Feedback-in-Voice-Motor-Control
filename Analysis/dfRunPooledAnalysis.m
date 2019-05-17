function dfRunPooledAnalysis()
% dfRunPooledAnalysis() opens run result files from multiple subjects and 
% multiple runs and pools the results for further statistics and plotting.
% To run this function, it is required that you have created a GenConfig 
% function within the results folder and generated a PooledConfig MATLAB 
% data file. PooledConfig is a specific configuration of the order of 
% subjects and runs to load. This keeps things neat so that 
% dfRunPooledAnalysis() can be used for more than one set of recordings.
%
% Different PooledConfig files can be selected by inputing the name of the
% pooled data set at line 20.
% 
% Dependent on the following packages from the 'MATLAB-Toolboxes' Repo
% -swtest
%
% Requires the Signal Processing Toolbox

close all
pA.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control'; 
pA.pAnalysis     = 'DRF_Aud'; % Change this name to load different pooled data sets Ex: SfN2017, LarynxPos

dirs               = dfDirs(pA.project);
dirs.SavResultsDir = fullfile(dirs.Results, 'Pooled Analyses', pA.pAnalysis);
dirs.PooledConfigF = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'PooledConfig.mat']);

% Can we find the Pooled Config File?
if exist(dirs.PooledConfigF, 'file') == 0
    fprintf('\nERROR: Pooled Config File %s does not exist! Please create it with a GenConfig Function\n', dirs.PooledConfigF)
    return
else
    % Load the configuration file. Should return a data structure 'cF'
    load(dirs.PooledConfigF)
end 

pA.participants  = cF.participants; % List of multiple participants.
pA.numPart       = length(pA.participants);
pA.runs          = cF.runs;         % All runs to consider 
pA.numRuns       = length(pA.runs);
pA.cond          = cF.cond;         % Conditions to test against
pA.numCond       = length(pA.cond); 
pA.condVar       = cF.condVar;      % Variable to test the condition
pA.testExt       = cF.testExt;
pA.pubCond       = cF.pubCond;

pA.pltNameMVi    = cell(pA.numPart, 1);
pA.pltNameMVm    = [pA.pAnalysis 'MeanSubj' pA.testExt];

% Load all saved results and order into a large data structure
allDataStr = []; % (numPart x numRuns)
for ii = 1:pA.numPart
    participant = pA.participants{ii};
    fprintf('Loading Runs for %s\n', participant)
    
    pA.pltNameMVi{ii} = [pA.pAnalysis participant pA.testExt];
    subjRes  = [];
    for jj = 1:pA.numRuns
        run              = pA.runs{jj};
        dirs.SavFileDir  = fullfile(dirs.Results, participant, run);                      % Where results are saved
        dirs.SavFile     = fullfile(dirs.SavFileDir, [participant run 'ResultsDRF.mat']); % Run Results file to load

        if exist(dirs.SavFile, 'file') == 0
            error('ERROR: The Results file for Run %s does not exist yet\n', run)
        else   
            load(dirs.SavFile)
            % Returns a results struture of 'res'
        end
        
        subjRes = cat(2, subjRes, res);  
    end
    allDataStr = cat(1, allDataStr, subjRes);
end

allSubjRes         = initSortedStruct(pA);
allSubjRes.subject = 'Mean Participant Response';
allSubjRes.curSess = allSubjRes.subject;
allSubjRes.cond    = pA.cond;
allSubjRes.pubCond = pA.pubCond;

[tossTrialTracker, tVN] = initTossedTrialTracker();

for ii = 1:pA.numPart
    participant = pA.participants{ii};
    fprintf('Sorting task conditions for %s\n', participant)
    
    sortStruc         = initSortedStruct(pA);
    sortStruc.subject = ['Participant ' num2str(ii)]; % Pooled Analysis Name
    
    sortStruc.curSess = sortStruc.subject;
    sortStruc.cond    = pA.cond;
    sortStruc.pubCond = pA.pubCond;
 
    for jj = 1:pA.numRuns
        curRes = allDataStr(ii, jj);

        sortStruc.studyID = curRes.subject; % Study ID
        sortStruc.expType = curRes.expType;
        sortStruc.gender  = curRes.gender;
        sortStruc.age     = round(curRes.age, 1);
        
        [tossCounts, tossTrialTracker, autoMiss] = combineTossedTrialTracker(curRes, tossTrialTracker);
        
        if ~isempty(autoMiss)
            curRes = adjustCurRes(curRes, autoMiss);
        end
        
        sortStruc  = combineCondTrials(pA, curRes, sortStruc, tossCounts);       
        allSubjRes = combineCondTrials(pA, curRes, allSubjRes, tossCounts);
    end
        
    sortStruc = meanCondTrials(pA, sortStruc);
    sortStruc.pltName = pA.pltNameMVi{ii};
    
    % Organize Stat Tables
    sortStruc  = catStatTableObs(pA, sortStruc, sortStruc);
    allSubjRes = catStatTableObs(pA, allSubjRes, sortStruc);
    sortStruc.statTable = packStatTable(sortStruc);
    
    pooledRunStr(ii)   = sortStruc;
    
    allSubjRes.expType    = sortStruc.expType;
    allSubjRes.gender{ii} = sortStruc.gender;
    allSubjRes.age(ii)    = sortStruc.age;
end

allSubjRes = catSubjMeans(pA, allSubjRes, pooledRunStr);
allSubjRes = meanCondTrials(pA, allSubjRes);
allSubjRes.pltName  = pA.pltNameMVm;

allSubjRes.statTable       = packStatTable(allSubjRes);
allSubjRes.statTableSingle = packStatTableSingle(pooledRunStr);

% Organize and Print the Stats of the Demographics included in this study
organizeAndPrintDemographicStats(dirs, allSubjRes);

% Organize and Save the Table of Excluded Trials
organizeAndSaveExcludedTrialTable(dirs, pA, allSubjRes, tossTrialTracker, tVN, 1)

% organizePooledResultsForFrank(dirs, allSubjRes)

% Save the Pooled Results
dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'ResultsDRF.mat']);
fprintf('Saving Pooled Analysis for %s\n', pA.pAnalysis)
save(dirs.SavResultsFile, 'pooledRunStr', 'allSubjRes')

%Apply Stats as appropriate
if strcmp(pA.pAnalysis, 'MaskingStudy')
    StatsOrg_MaskingNoiseStudy(dirs, pA, allSubjRes);
elseif strcmp(pA.pAnalysis, 'DRF_Som')
    StatsOrg_DRF_Som(dirs, pA, allSubjRes);
    timeSeriesDiffAnalysis(dirs, pA, allSubjRes)
    drawSubjRespVarDists(dirs, pooledRunStr)
elseif strcmp(pA.pAnalysis, 'DRF_Aud')
    drawSubjRespVarDists(dirs, pooledRunStr)
    StatsOrg_DRF_Aud(dirs, pA, allSubjRes);
end
end

function sortStr = initSortedStruct(pA)
% sortStr = initSortedStruct(numCond) initializes the structure that will
% store the pooled results for each subject, or group of subjects. It is
% created to have different sizes, based on the number of conditions that
% are being tested against. I think this should generalize to subconditions
% of conditions, or two condition crossing, but I have not tested that, and
% currently the above scripts only consider one condition to test against. 

pAnalysis = pA.pAnalysis;
numCond = pA.numCond;

% Basic info about the session, the recordings, the subjects
sortStr.pAnalysis = pAnalysis;
sortStr.expType = [];
sortStr.subject = [];
sortStr.gender  = [];
sortStr.f0      = [];
sortStr.age     = [];
sortStr.curSess = [];
sortStr.studyID = [];
sortStr.runs    = cell(numCond, 1);
sortStr.runf0b  = cell(numCond, 1);
sortStr.f0b     = zeros(numCond, 1);

sortStr.AudFB   = cell(numCond, 1);

sortStr.allContTrials    = [];
sortStr.numContTrialsFin = 0;
sortStr.allPertTrials    = cell(numCond, 1);
sortStr.numPertTrialsFin = zeros(numCond, 1);

sortStr.secTime         = [];
sortStr.audioMf0SecPert = cell(numCond, 1);
sortStr.audioMf0SecCont = [];
sortStr.audioHf0SecPert = cell(numCond, 1);
sortStr.audioHf0SecCont = [];

sortStr.secTimeP        = [];
sortStr.sensorPSec      = [];

sortStr.audioMf0MeanPert = cell(numCond, 1);
sortStr.audioMf0MeanCont = [];
sortStr.audioHf0MeanPert = cell(numCond, 1);
sortStr.audioHf0MeanCont = [];
sortStr.sensorPMean      = [];

sortStr.tossedAll        = 0;
sortStr.tossedLate       = 0;
sortStr.tossedBreak      = 0;
sortStr.tossedMisCalc    = 0;
sortStr.tossedManual     = 0;
sortStr.tossedAutoMiss   = 0;

sortStr.respVarSingle    = cell(numCond, 1);
sortStr.respVarM         = cell(numCond, 1);

sortStr.obvSubj          = {};
sortStr.obvAge           = [];
sortStr.obvGender        = {};
sortStr.obvAudFB         = {};
sortStr.obvRespVar       = [];
end

function [tossTrialTracker, tVN] = initTossedTrialTracker()

tossTrialTracker              = []; %structure of keeping of which trials and when
tossTrialTracker.curSess      = {};
tossTrialTracker.tossedTrials = {};
tossTrialTracker.manuallyExcl = {};
tossTrialTracker.autoMiss     = {};
tossTrialTracker.perCaught    = {};

% Tossed Trials Table Variable Names
tVN = {'CurSess', 'AutoExcluded', 'ManuallyExcluded', 'MissedByAuto', 'PercentCaught'};
end

function polRes = combineCondTrials(pA, curRes, polRes, tossT)
AD = curRes.audioDynamics; % Audio Dynamics
PD = curRes.presSDsv;      % Pressure Dynamics

whichCondAr = strcmp(pA.cond, eval(pA.condVar));
wC          = find(whichCondAr == 1);            % Which Condition?

polRes.runs{wC}   = cat(1, polRes.runs{wC}, {curRes.run});
polRes.runf0b{wC} = cat(1, polRes.runf0b{wC}, curRes.f0b);

polRes.AudFB{wC}  = cat(1, polRes.AudFB{wC}, {curRes.AudFB});

polRes.allContTrials     = cat(1, polRes.allContTrials, curRes.numContTrialsFin);
polRes.allPertTrials{wC} = cat(1, polRes.allPertTrials{wC}, curRes.numPertTrialsFin);

polRes.secTime             = curRes.secTime;
polRes.audioMf0SecPert{wC} = cat(2, polRes.audioMf0SecPert{wC}, curRes.audioMf0SecPert);
polRes.audioHf0SecPert{wC} = cat(2, polRes.audioHf0SecPert{wC}, curRes.audioHf0SecPert);
polRes.audioMf0SecCont     = cat(2, polRes.audioMf0SecCont, curRes.audioMf0SecCont);
polRes.audioHf0SecCont     = cat(2, polRes.audioHf0SecCont, curRes.audioHf0SecCont);

polRes.secTimeP            = PD.timeSec;
polRes.sensorPSec          = cat(2, polRes.sensorPSec, PD.sensorSec);

polRes.tossedAll      = polRes.tossedAll + tossT.A;       % Total Automatic Excluded Trials
polRes.tossedLate     = polRes.tossedLate + tossT.L;      % Late Start
polRes.tossedBreak    = polRes.tossedBreak + tossT.B;     % Voice Break
polRes.tossedMisCalc  = polRes.tossedMisCalc + tossT.C;   % f0 Miscalc
polRes.tossedManual   = polRes.tossedManual + tossT.M;    % Total Manual Excluded Trials
polRes.tossedAutoMiss = polRes.tossedAutoMiss + tossT.aM; % Trials Manually removed, but missed by auto methods.
end

function [tossCounts, tossTrialTracker, autoMiss] = combineTossedTrialTracker(curRes, tossTrialTracker)

% AUTO: Identify trials that were excluded by automated methods 
[tossCounts, aTossTrials, aTossDetails] = identifyTossedTrials(curRes);

% MANUAL: Identify trials that were excluded by manual methods
[mTossTrials, mTossDetails] = unpackManuallyExcludedTrials(curRes.incTrialInfo);

% COMPARE AUTO and MANUAL
[autoMiss, perCaught] = compareTossedTrials(aTossTrials, mTossTrials);
numManual = length(mTossTrials);
numMissed = length(autoMiss);
tossCounts.M  = numManual;
tossCounts.aM = numMissed;

tossTrialTracker.curSess      = cat(1, tossTrialTracker.curSess, curRes.curSess);
tossTrialTracker.tossedTrials = cat(1, tossTrialTracker.tossedTrials, {aTossDetails});
tossTrialTracker.manuallyExcl = cat(1, tossTrialTracker.manuallyExcl, {mTossDetails});
tossTrialTracker.autoMiss     = cat(1, tossTrialTracker.autoMiss, {num2str(autoMiss)});
tossTrialTracker.perCaught    = cat(1, tossTrialTracker.perCaught, {perCaught}); 
end

function [tossCounts, aTossTrials, aTossDetails] = identifyTossedTrials(curRes)

tT = curRes.removedTrialTracker;
aTossDetails = [];
if ~isempty(tT)
    [tossAll, ~] = size(tT);
    
    tossLate  = strcmp(tT(:,2), 'Participant started too late!!');
    tossBreak = strcmp(tT(:,2), 'Participant had a voice break!!');
    tossCalc  = strcmp(tT(:,2), 'Miscalculated pitch Trace');

    tossCounts.A = sum(tossAll); 
    tossCounts.L = sum(tossLate);
    tossCounts.B = sum(tossBreak);
    tossCounts.C = sum(tossCalc);
    
    [aTossDetails, LSTrialNums] = writeTossedDetails(tT, tossLate, 'LS', aTossDetails);
    [aTossDetails, VBTrialNums] = writeTossedDetails(tT, tossBreak, 'VB', aTossDetails);
    [aTossDetails, MCTrialNums] = writeTossedDetails(tT, tossCalc, 'MC', aTossDetails);
    aTossTrials  = [LSTrialNums VBTrialNums MCTrialNums];
else
    tossCounts.A = 0;
    tossCounts.L = 0;
    tossCounts.B = 0;
    tossCounts.C = 0;
    aTossDetails = '';
    aTossTrials  = [];
end

end

function [tDetails, trialNums] = writeTossedDetails(tT, logicToss, note, tDetails)
% [tDetails, trialNums] = writeTossedDetails(tT, logicToss, note, tDetails)
% records which trials were lost and for reason for use in a table later on
% which details all trials considered.

numTrials = sum(logicToss);
curTrials = tT(logicToss);

trialNums = [];
for i = 1:numTrials
    curTrial = curTrials{i};
    spc = find(curTrial == ' ');
    curTrialNum = curTrial(spc+1:end);
    
    tDetails  = cat(2, tDetails, [curTrialNum '(' note ') ']);
    trialNums = cat(2, trialNums, str2double(curTrialNum));
end
end

function [mTossTrials, mTossDetails] = unpackManuallyExcludedTrials(incTrialInfo)
% [mTossTrials, mTossDetails] = unpackManuallyExcludedTrials(incTrialInfo)
% opens the notes we made regarding which trials to manually exclude. These
% notes are later compared against the automated methods
mTossTrials  = [];
mTossDetails = [];

if ~isempty(incTrialInfo)
    excludedTrials     = find([incTrialInfo.include] == 0);
    excludedTrialNotes = [incTrialInfo(excludedTrials).trialNote];
    numExcludedTrials  = length(excludedTrials);

    for ii = 1:numExcludedTrials
        curExcTrial = excludedTrials(ii);

        mTossTrials  = cat(2, mTossTrials, curExcTrial);
        mTossDetails = cat(2, mTossDetails, [num2str(curExcTrial) ' (' excludedTrialNotes{ii} ') ']);
    end
end
end

function [autoMiss, perCaughtStr] = compareTossedTrials(aTossTrials, mTossTrials)

numM = length(mTossTrials);
if numM == 0
    perCaughtStr = 'N/A';
else
    autoCorrect  = intersect(aTossTrials, mTossTrials);
    numCorrect   = length(autoCorrect);
    perCaught    = round(100*(numCorrect/numM));
    perCaughtStr = [num2str(perCaught) '%'];
end

autoMiss = setdiff(mTossTrials, aTossTrials);
end

function curRes = adjustCurRes(curRes, autoMiss)

keptElementsCont = curRes.allIdxFin(curRes.contIdxFin);
keptElementsPert = curRes.allIdxFin(curRes.pertIdxFin);

[~, contIdx2Remove] = intersect(keptElementsCont, autoMiss);
[~, pertIdx2Remove] = intersect(keptElementsPert, autoMiss);

%Re-Adjust curRes for those extra trials we threw away.
curRes.audioMf0SecCont(:, contIdx2Remove,:) = [];
curRes.audioMf0SecPert(:, pertIdx2Remove,:) = [];
curRes.presSDsv.sensorSec(:, pertIdx2Remove,:)      = [];

[~, curRes.numContTrialsFin, ~] = size(curRes.audioMf0SecCont);
[~, curRes.numPertTrialsFin, ~] = size(curRes.audioMf0SecPert);
end

function polRes = meanCondTrials(pA, polRes)

polRes.numContTrialsFin = sum(polRes.allContTrials);
polRes.audioMf0MeanCont = meanSecData(polRes.audioMf0SecCont);
polRes.audioHf0MeanCont = meanSecData(polRes.audioHf0SecCont);

if strcmp(polRes.expType, 'Somatosensory Perturbation_Perceptual')
    polRes.sensorPMean = meanSecData(polRes.sensorPSec);
else
    polRes.sensorPMean = zeros(size(polRes.audioMf0SecPert{1}));
end

for wC = 1:pA.numCond
    polRes.f0b(wC)              = mean(polRes.runf0b{wC});
    
    polRes.numPertTrialsFin(wC) = sum(polRes.allPertTrials{wC});
    polRes.audioMf0MeanPert{wC} = meanSecData(polRes.audioMf0SecPert{wC});
    polRes.audioHf0MeanPert{wC} = meanSecData(polRes.audioHf0SecPert{wC});
    
    if strcmp(polRes.expType, 'Somatosensory Perturbation_Perceptual')
        audioDynamicsAllTraces = InflationResponseSingle(polRes.secTime, polRes.audioMf0SecPert{wC});
        audioDynamicsMeanTrace = InflationResponse(polRes.secTime, polRes.audioMf0MeanPert{wC});
    elseif strcmp(polRes.expType, 'Auditory Perturbation_Perceptual')
        audioDynamicsAllTraces = PitchShiftReflexResponseSingle(polRes.secTime, polRes.audioMf0SecPert{wC});
        audioDynamicsMeanTrace = PitchShiftReflexResponse(polRes.secTime, polRes.audioMf0MeanPert{wC});
    end
    polRes.respVarSingle{wC} = audioDynamicsAllTraces;
    polRes.respVarM{wC} = audioDynamicsMeanTrace.respVarM;
end

% Identify Limits of the newly meaned data
lims = identifyLimits(polRes);
polRes.limitsAmean = lims.audioMean;
polRes.limitsMHmean = lims.audioMHMean;
polRes.limitsPmean = lims.presMean;

% Scale Pressure traces against the f0 traces
if strcmp(polRes.expType, 'Somatosensory Perturbation_Perceptual')
    [sensorPAdjust, InflDeflT] = scalePressureVals(polRes);
else
    sensorPAdjust = polRes.sensorPMean;
    InflDeflT     = [0.04 0.24 1.04 1.24];
end
polRes.sensorPAdjust = sensorPAdjust;
polRes.InflDeflT     = InflDeflT;

% Identify success of automatic trial exclusion methods
if polRes.tossedManual > 0
    perc = round(100*(1-(polRes.tossedAutoMiss/polRes.tossedManual)));
    polRes.autoSuccessPerc = num2str(perc);
else
    polRes.autoSuccessPerc = 'N/A';
end
end

function polRes = catStatTableObs(pA, polRes, sortStruc)

for wC = 1:pA.numCond
    polRes.obvSubj    = cat(1, polRes.obvSubj, sortStruc.studyID);
    polRes.obvAge     = cat(1, polRes.obvAge, sortStruc.age);
    polRes.obvGender  = cat(1, polRes.obvGender, sortStruc.gender);
    polRes.f0         = cat(1, polRes.f0, sortStruc.f0b(wC));
    polRes.obvAudFB   = cat(1, polRes.obvAudFB, sortStruc.AudFB{wC}{1});
    polRes.obvRespVar = cat(1, polRes.obvRespVar, round(sortStruc.respVarM{wC},2));
end
end

function allSubjRes = catSubjMeans(pA, allSubjRes, pooledRunStr)

numPartici = length(pooledRunStr);

audioMf0SecContAllSubj = [];
audioMf0SecPertAllSubj = cell(pA.numCond,1);
audioHf0SecContAllSubj = [];
audioHf0SecPertAllSubj = cell(pA.numCond,1);
for nP = 1:numPartici
    curSubj = pooledRunStr(nP);
    
    meanOnsetCont   = curSubj.audioMf0MeanCont(:,1);
    meanOffsetCont  = curSubj.audioMf0MeanCont(:,3);
    secCont(:,:,1) = meanOnsetCont;
    secCont(:,:,2) = meanOffsetCont;
    audioMf0SecContAllSubj = cat(2, audioMf0SecContAllSubj, secCont);
    
    meanOnsetCont   = curSubj.audioHf0MeanCont(:,1);
    meanOffsetCont  = curSubj.audioHf0MeanCont(:,3);
    secCont(:,:,1) = meanOnsetCont;
    secCont(:,:,2) = meanOffsetCont;
    audioHf0SecContAllSubj = cat(2, audioHf0SecContAllSubj, secCont);
    
    for wC = 1:pA.numCond
        meanOnsetPert   = curSubj.audioMf0MeanPert{wC}(:,1);
        meanOffsetPert  = curSubj.audioMf0MeanPert{wC}(:,3);
        secPert(:,:,1) = meanOnsetPert;
        secPert(:,:,2) = meanOffsetPert;
        audioMf0SecPertAllSubj{wC} = cat(2, audioMf0SecPertAllSubj{wC}, secPert);
        
        meanOnsetPert   = curSubj.audioHf0MeanPert{wC}(:,1);
        meanOffsetPert  = curSubj.audioHf0MeanPert{wC}(:,3);
        secPert(:,:,1) = meanOnsetPert;
        secPert(:,:,2) = meanOffsetPert;
        audioHf0SecPertAllSubj{wC} = cat(2, audioHf0SecPertAllSubj{wC}, secPert);
    end
end

allSubjRes.audioMf0SecCont = audioMf0SecContAllSubj;
allSubjRes.audioMf0SecPert = audioMf0SecPertAllSubj;
allSubjRes.audioHf0SecCont = audioHf0SecContAllSubj;
allSubjRes.audioHf0SecPert = audioHf0SecPertAllSubj;
end

function meanData = meanSecData(secData)

OnsetSecs  = secData(:,:,1);
OffsetSecs = secData(:,:,2);
[~, numTrial] = size(OnsetSecs);

meanOnset  = nanmean(OnsetSecs, 2);
meanOffset = nanmean(OffsetSecs, 2);

stdOnset   = nanstd(OnsetSecs, 0, 2);
stdOffset  = nanstd(OffsetSecs, 0, 2);

SEMOnset   = stdOnset/sqrt(numTrial);  % Standard Error
SEMOffset  = stdOffset/sqrt(numTrial); % Standard Error

% NCIOnset   = 1.96*SEMOnset;  % 95% Confidence Interval
% NCIOffset  = 1.96*SEMOffset; % 95% Confidence Interval

meanData   = [meanOnset SEMOnset meanOffset SEMOffset];
end

function audioDynamics_Somato = InflationResponse(secTime, secAudioMean)
% [respVar, respVarm, respVarSD, InflaStimVar] = InflationResponse(secTime, secAudio)
% Identifies the relevant pitch contour characteristics that are important
% for deciding how a participant responded to the inflation of the balloon
% during production. iR is a structure representing the result variables
% from studying the inflation response (iR). The prefix letter denotes
% whether the variable is a index (i), a time (t), or a value (v). 
%
% secTime:  Vector of time points corresponding to the sectioned data (numSamp)
% secAudio: 3D mat of sectioned audio (numSamp x numTrial x event)
%           The 1st 3D layer are Onset Sections
%           The 2nd 3D later are Offset Sections
%
% respVar: Matrix of per trial iR results (numTrial x 4)
%          respVar(:,1) = Time of the minimum f0 value in the sec
%          respVar(:,2) = Minimum f0 value in sec (stim magnitude)
%          respVar(:,3) = Value of f0 at end of sec (response magnitude)
%          respVar(:,4) = ABS percent change of stim and response
%          magnitude (response percentage)
% respVarM:    Vector of mean trial values from respVarm (1x4)
% respVarSD:   Vector of standard deviation of the trial values from respVar (1x4)
% InflaSimVar: Values of the mean time at stim magnitude and the mean stim magnitude

[numSamp, ~] = size(secAudioMean); % Size of the data we are dealing with

ir = initInflationResponseStruct(); % Initialize the structure that handles the variable calculations
ir.time     = secTime;              % Time Interval for the sectioned trials (-0.5->1.0s)
ir.onset    = secAudioMean(:, 1);   % f0 Trace sectioned around pert Onset.

ir.iAtOnset = find(ir.time == 0);
ir.tAtOnset = 0;                     % duh
ir.vAtOnset = ir.onset(ir.iAtOnset); % f0 value at t = 0

ir.iPostOnsetR = find(0 < ir.time & .20 >= ir.time); % Range of indices between t > 0ms and t =< 200ms;
[minOn, minIdx] = min(ir.onset(ir.iPostOnsetR));     % Minimum f0 val within PostOnsetR

% StimMag
ir.iAtMin  = ir.iPostOnsetR(minIdx);       % Indice of the min f0 value following trigger
ir.tAtMin  = ir.time(ir.iAtMin);           % Time at min f0 value following trigger
ir.vAtMin  = minOn;                        % Min f0 value in PostOnsetR
ir.stimMag = abs(ir.vAtMin - ir.vAtOnset); % Distance traveled from onset to min value

% RespMag
ir.iAtResp = numSamp;                % Last index in section
ir.tAtResp = ir.time(ir.iAtResp);    % Time Value when participant 'fully responded' (1.0s)
ir.vAtResp = ir.onset(ir.iAtResp);   % f0 value when participant 'fully responded'
ir.respMag = ir.vAtResp - ir.vAtMin; % Distance traveled from min f0 value to response f0 value

% RespPer
ir.respPer = 100*(ir.respMag/ir.stimMag); % Percent change from stimMag to respMag 

% Add to the audioDynamics struct
respVarM = [ir.tAtMin ir.stimMag ir.respMag ir.respPer];
audioDynamics_Somato.respVarM = respVarM;
end

function audioDynamics_Somato = InflationResponseSingle(secTime, secAudioMean)
% [respVar, respVarm, respVarSD, InflaStimVar] = InflationResponse(secTime, secAudio)
% Identifies the relevant pitch contour characteristics that are important
% for deciding how a participant responded to the inflation of the balloon
% during production. iR is a structure representing the result variables
% from studying the inflation response (iR). The prefix letter denotes
% whether the variable is a index (i), a time (t), or a value (v). 
%
% secTime:  Vector of time points corresponding to the sectioned data (numSamp)
% secAudio: 3D mat of sectioned audio (numSamp x numTrial x event)
%           The 1st 3D layer are Onset Sections
%           The 2nd 3D later are Offset Sections
%
% respVar: Matrix of per trial iR results (numTrial x 4)
%          respVar(:,1) = Time of the minimum f0 value in the sec
%          respVar(:,2) = Minimum f0 value in sec (stim magnitude)
%          respVar(:,3) = Value of f0 at end of sec (response magnitude)
%          respVar(:,4) = ABS percent change of stim and response
%          magnitude (response percentage)
% respVarM:    Vector of mean trial values from respVarm (1x4)
% respVarSD:   Vector of standard deviation of the trial values from respVar (1x4)
% InflaSimVar: Values of the mean time at stim magnitude and the mean stim magnitude

[numSamp, numTrial, ~] = size(secAudioMean); % Size of the data we are dealing with

audioDynamics_Somato.tAtMin  = [];
audioDynamics_Somato.stimMag = [];
audioDynamics_Somato.respMag = [];
audioDynamics_Somato.respPer = [];
for ii = 1:numTrial
    ir = initInflationResponseStruct(); % Initialize the structure that handles the variable calculations
    ir.time     = secTime;              % Time Interval for the sectioned trials (-0.5->1.0s)
    ir.onset    = secAudioMean(:, ii, 1);   % f0 Trace sectioned around pert Onset.

    ir.iAtOnset = find(ir.time == 0);
    ir.tAtOnset = 0;                     % duh
    ir.vAtOnset = ir.onset(ir.iAtOnset); % f0 value at t = 0

    ir.iPostOnsetR = find(0 < ir.time & .20 >= ir.time); % Range of indices between t > 0ms and t =< 200ms;
    [minOn, minIdx] = min(ir.onset(ir.iPostOnsetR));     % Minimum f0 val within PostOnsetR

    % StimMag
    ir.iAtMin  = ir.iPostOnsetR(minIdx);       % Indice of the min f0 value following trigger
    ir.tAtMin  = ir.time(ir.iAtMin);           % Time at min f0 value following trigger
    ir.vAtMin  = minOn;                        % Min f0 value in PostOnsetR
    ir.stimMag = abs(ir.vAtMin - ir.vAtOnset); % Distance traveled from onset to min value

    % RespMag
    ir.iAtResp = numSamp;                % Last index in section
    ir.tAtResp = ir.time(ir.iAtResp);    % Time Value when participant 'fully responded' (1.0s)
    ir.vAtResp = ir.onset(ir.iAtResp);   % f0 value when participant 'fully responded'
    ir.respMag = ir.vAtResp - ir.vAtMin; % Distance traveled from min f0 value to response f0 value

    % RespPer
    ir.respPer = 100*(ir.respMag/ir.stimMag); % Percent change from stimMag to respMag
    
    audioDynamics_Somato.tAtMin  = cat(1, audioDynamics_Somato.tAtMin, ir.tAtMin);
    audioDynamics_Somato.stimMag = cat(1, audioDynamics_Somato.stimMag, ir.stimMag);
    audioDynamics_Somato.respMag = cat(1, audioDynamics_Somato.respMag, ir.respMag);
    audioDynamics_Somato.respPer = cat(1, audioDynamics_Somato.respPer, ir.respPer);
end

audioDynamics_Somato.tAtMinM   = round(mean(audioDynamics_Somato.tAtMin),3);
audioDynamics_Somato.tAtMinSD  = round(std(audioDynamics_Somato.tAtMin),3);
audioDynamics_Somato.stimMagM  = round(mean(audioDynamics_Somato.stimMag),2);
audioDynamics_Somato.stimMagSD = round(std(audioDynamics_Somato.stimMag),2);
audioDynamics_Somato.respMagM  = round(mean(audioDynamics_Somato.respMag),2);
audioDynamics_Somato.respMagSD = round(std(audioDynamics_Somato.respMag),2);
audioDynamics_Somato.respPerM  = round(mean(audioDynamics_Somato.respPer),2);
audioDynamics_Somato.respPerSD = round(std(audioDynamics_Somato.respPer),2);
end

function ir = initInflationResponseStruct()

ir.numSamp  = [];
ir.numTrial = [];
ir.time     = [];
ir.onset    = [];

ir.iAtOnset = []; % Index where t = 0
ir.tAtOnset = []; % Time at t = 0
ir.vAtOnset = []; % f0 value at t = 0

ir.iPostOnsetR = []; % Range of indices between t = 0ms and t = 200ms;
ir.iAtMin      = []; % Index at min f0 value in PostOnsetR
ir.tAtMin      = []; % Time at min f0 value in PostOnsetR
ir.vAtMin      = []; % Min f0 value in PostOnsetR
ir.stimMag     = []; % ir.vAtMin - ir.vAtOnset ..in a perfect world vAtOnset = 0

ir.iAtResp = []; % Index of f0 value when participant 'fully' responded...right now = last value in section
ir.tAtResp = []; % Time at f0 value when participant 'fully' responded
ir.vAtResp = []; % f0 value when participant 'fully' responded 
ir.respMag = []; % vAtResp - vAtMin   ...distance traveled
ir.respPer = []; % Percent change from stimMag to respMag
end

function audioDynamics_Audio = PitchShiftReflexResponse(secTime, secAudioMean)
% [respVar, respVarm, respVarSD, InflaStimVar] = InflationResponse(secTime, secAudio)
% Identifies the relevant pitch contour characteristics that are important
% for deciding how a participant responded to the inflation of the balloon
% during production. iR is a structure representing the result variables
% from studying the inflation response (iR). The prefix letter denotes
% whether the variable is a index (i), a time (t), or a value (v). 
%
% secTime:  Vector of time points corresponding to the sectioned data (numSamp)
% secAudio: 3D mat of sectioned audio (numSamp x numTrial x event)
%           The 1st 3D layer are Onset Sections
%           The 2nd 3D later are Offset Sections
%
% respVar: Matrix of per trial iR results (numTrial x 4)
%          respVar(:,1) = Time of the minimum f0 value in the sec
%          respVar(:,2) = Minimum f0 value in sec (stim magnitude)
%          respVar(:,3) = Value of f0 at end of sec (response magnitude)
%          respVar(:,4) = ABS percent change of stim and response
%          magnitude (response percentage)
% respVarM:    Vector of mean trial values from respVarm (1x4)
% respVarSD:   Vector of standard deviation of the trial values from respVar (1x4)
% InflaSimVar: Values of the mean time at stim magnitude and the mean stim magnitude

[numSamp, ~] = size(secAudioMean); % Size of the data we are dealing with

ir = initInflationResponseStruct(); % Initialize the structure that handles the variable calculations
ir.time     = secTime;              % Time Interval for the sectioned trials (-0.5->1.0s)
ir.onset    = secAudioMean(:, 1);   % f0 Trace sectioned around pert Onset.

ir.iAtOnset = find(ir.time == 0);
ir.tAtOnset = 0;                     % duh
ir.vAtOnset = ir.onset(ir.iAtOnset); % f0 value at t = 0

ir.iPostOnsetR = find(0 < ir.time & .20 >= ir.time); % Range of indices between t > 0ms and t =< 200ms;
[minOn, minIdx] = min(ir.onset(ir.iPostOnsetR));     % Minimum f0 val within PostOnsetR

% StimMag
ir.iAtMin  = ir.iPostOnsetR(minIdx);       % Indice of the min f0 value following trigger
ir.tAtMin  = ir.time(ir.iAtMin);           % Time at min f0 value following trigger
ir.vAtMin  = minOn;                        % Min f0 value in PostOnsetR
ir.stimMag = abs(-100);                    % Distance traveled from onset to min value (default is 100 cents)

% RespMag
ir.iAtResp = numSamp;                % Last index in section
ir.tAtResp = ir.time(ir.iAtResp);    % Time Value when participant 'fully responded' (1.0s)
ir.vAtResp = ir.onset(ir.iAtResp);   % f0 value when participant 'fully responded'
ir.respMag = ir.vAtResp - ir.vAtMin; % Distance traveled from min f0 value to response f0 value

% RespPer
ir.respPer = 100*(ir.respMag/ir.stimMag); % Percent change from stimMag to respMag 

% Add to the audioDynamics struct
respVarM = [ir.tAtMin ir.stimMag ir.respMag ir.respPer];
audioDynamics_Audio.respVarM = respVarM;
end

function audioDynamics_Audio = PitchShiftReflexResponseSingle(secTime, secAudioMean)
% [respVar, respVarm, respVarSD, InflaStimVar] = InflationResponse(secTime, secAudio)
% Identifies the relevant pitch contour characteristics that are important
% for deciding how a participant responded to the inflation of the balloon
% during production. iR is a structure representing the result variables
% from studying the inflation response (iR). The prefix letter denotes
% whether the variable is a index (i), a time (t), or a value (v). 
%
% secTime:  Vector of time points corresponding to the sectioned data (numSamp)
% secAudio: 3D mat of sectioned audio (numSamp x numTrial x event)
%           The 1st 3D layer are Onset Sections
%           The 2nd 3D later are Offset Sections
%
% respVar: Matrix of per trial iR results (numTrial x 4)
%          respVar(:,1) = Time of the minimum f0 value in the sec
%          respVar(:,2) = Minimum f0 value in sec (stim magnitude)
%          respVar(:,3) = Value of f0 at end of sec (response magnitude)
%          respVar(:,4) = ABS percent change of stim and response
%          magnitude (response percentage)
% respVarM:    Vector of mean trial values from respVarm (1x4)
% respVarSD:   Vector of standard deviation of the trial values from respVar (1x4)
% InflaSimVar: Values of the mean time at stim magnitude and the mean stim magnitude

[numSamp, numTrial, ~] = size(secAudioMean); % Size of the data we are dealing with

audioDynamics_Audio.stimMag = [];
audioDynamics_Audio.respMag = [];
audioDynamics_Audio.respPer = [];
for ii = 1:numTrial
    ir = initInflationResponseStruct(); % Initialize the structure that handles the variable calculations
    ir.time     = secTime;              % Time Interval for the sectioned trials (-0.5->1.0s)
    ir.onset    = secAudioMean(:, ii, 1);   % f0 Trace sectioned around pert Onset.

    ir.iAtOnset = find(ir.time == 0);
    ir.tAtOnset = 0;                     % duh
    ir.vAtOnset = ir.onset(ir.iAtOnset); % f0 value at t = 0

    ir.iPostOnsetR = find(0 < ir.time & .20 >= ir.time); % Range of indices between t > 0ms and t =< 200ms;
    [minOn, minIdx] = min(ir.onset(ir.iPostOnsetR));     % Minimum f0 val within PostOnsetR

    % StimMag
    ir.iAtMin  = ir.iPostOnsetR(minIdx);       % Indice of the min f0 value following trigger
    ir.tAtMin  = ir.time(ir.iAtMin);           % Time at min f0 value following trigger
    ir.vAtMin  = minOn;                        % Min f0 value in PostOnsetR
    ir.stimMag = abs(-100);                    % Distance traveled from onset to min value (default is 100 cents)

    % RespMag
    ir.iAtResp = numSamp;                % Last index in section
    ir.tAtResp = ir.time(ir.iAtResp);    % Time Value when participant 'fully responded' (1.0s)
    ir.vAtResp = ir.onset(ir.iAtResp);   % f0 value when participant 'fully responded'
    ir.respMag = ir.vAtResp - ir.vAtMin; % Distance traveled from min f0 value to response f0 value

    % RespPer
    ir.respPer = 100*(ir.respMag/ir.stimMag); % Percent change from stimMag to respMag
    
    audioDynamics_Audio.stimMag = cat(1, audioDynamics_Audio.stimMag, ir.stimMag);
    audioDynamics_Audio.respMag = cat(1, audioDynamics_Audio.respMag, ir.respMag);
    audioDynamics_Audio.respPer = cat(1, audioDynamics_Audio.respPer, ir.respPer);
end

audioDynamics_Audio.stimMagM  = round(mean(audioDynamics_Audio.stimMag),2);
audioDynamics_Audio.stimMagSD = round(std(audioDynamics_Audio.stimMag),2);
audioDynamics_Audio.respMagM  = round(mean(audioDynamics_Audio.respMag),2);
audioDynamics_Audio.respMagSD = round(std(audioDynamics_Audio.respMag),2);
audioDynamics_Audio.respPerM  = round(mean(audioDynamics_Audio.respPer),2);
audioDynamics_Audio.respPerSD = round(std(audioDynamics_Audio.respPer),2);
end

function lims = identifyLimits(ss)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f0 Traces Limits
mf0MeanPert = ss.audioMf0MeanPert;
hf0MeanPert = ss.audioHf0MeanPert;
numCond     = length(mf0MeanPert);

setMupBoundSec = zeros(numCond, 1);
setMlwBoundSec = zeros(numCond, 1);
setHupBoundSec = zeros(numCond, 1);
setHlwBoundSec = zeros(numCond, 1);
for ii = 1:numCond
    audioMMean = mf0MeanPert{ii};    
    audioHMean = hf0MeanPert{ii};    

    [~, Imax] = max(audioMMean(:,1)); %Max Pert Onset
    upBoundOn = round(audioMMean(Imax,1) + audioMMean(Imax,2) + 10);
    [~, Imin] = min(audioMMean(:,1)); %Min Pert Onset
    lwBoundOn = round(audioMMean(Imin,1) - audioMMean(Imin,2) - 10);

    [~, Imax] = max(audioMMean(:,3)); %Max Pert Offset
    upBoundOf = round(audioMMean(Imax,3) + audioMMean(Imax,4) + 10);
    [~, Imin] = min(audioMMean(:,3)); %Min Pert Offset
    lwBoundOf = round(audioMMean(Imin,3) - audioMMean(Imin,4) - 10);

    if upBoundOn > upBoundOf
        upBoundSec = upBoundOn;
    else
        upBoundSec = upBoundOf;
    end

    if lwBoundOn < lwBoundOf
        lwBoundSec = lwBoundOn;
    else
        lwBoundSec = lwBoundOf;
    end
    
    setMupBoundSec(ii) = upBoundSec;
    setMlwBoundSec(ii) = lwBoundSec;
    
    [~, Imax] = max(audioHMean(:,1)); %Max Pert Onset
    upBoundOn = round(audioHMean(Imax,1) + audioHMean(Imax,2) + 10);
    [~, Imin] = min(audioHMean(:,1)); %Min Pert Onset
    lwBoundOn = round(audioHMean(Imin,1) - audioHMean(Imin,2) - 10);

    [~, Imax] = max(audioHMean(:,3)); %Max Pert Offset
    upBoundOf = round(audioHMean(Imax,3) + audioHMean(Imax,4) + 10);
    [~, Imin] = min(audioHMean(:,3)); %Min Pert Offset
    lwBoundOf = round(audioHMean(Imin,3) - audioHMean(Imin,4) - 10);

    if upBoundOn > upBoundOf
        upBoundSec = upBoundOn;
    else
        upBoundSec = upBoundOf;
    end

    if lwBoundOn < lwBoundOf
        lwBoundSec = lwBoundOn;
    else
        lwBoundSec = lwBoundOf;
    end
    
    setHupBoundSec(ii) = upBoundSec;
    setHlwBoundSec(ii) = lwBoundSec;
end

maxMUpBound = max(setMupBoundSec); % Max f0 Bound
minMLwBound = min(setMlwBoundSec); % Min f0 Bound 

maxHUpBound = max(setHupBoundSec); % Max f0 Bound
minHLwBound = min(setHlwBoundSec); % Min f0 Bound 

lims.audioMean = [-0.5 1.0 minMLwBound maxMUpBound];

maxMHUpBound = max([maxMUpBound maxHUpBound]);
minMHLwBound = min([minMLwBound minHLwBound]);

lims.audioMHMean = [-0.5 1.0 minMHLwBound maxMHUpBound];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pressure Traces Limits
maxPres = max(ss.sensorPMean(:,1)) + 0.5;
minPres = min(ss.sensorPMean(:,1)) - 0.1;

lims.presMean = [-0.5 1.0 minPres maxPres];
end

function [adjustPres, InflDeflT] = scalePressureVals(ss)

time        = ss.secTimeP;
presOnsetM  = ss.sensorPMean(:,1);
presOnsetE  = ss.sensorPMean(:,2);
presOffsetM = ss.sensorPMean(:,3);
presOffsetE = ss.sensorPMean(:,4);

limsA    = ss.limitsAmean;
limsAMin = limsA(3); limsAMax = limsA(4);
limsP    = ss.limitsPmean;
limsPMin = limsP(3); limsPMax = limsP(4);

% Scale pressure values against the 
m = (limsAMax - limsAMin)/(limsPMax - limsPMin);
b = limsAMin - m*limsPMin;

adjPresOnsetM  = (m*presOnsetM + b);
adjPresOnsetE  = (m*presOnsetE + b);
adjPresOffsetM = (m*presOffsetM + b);
adjPresOffsetE = (m*presOffsetE + b);

adjPresOnsetMin = min(adjPresOnsetM);
[~, maxInd]     = max(adjPresOnsetM);
[minInds]       = find(adjPresOnsetM > 0.98*adjPresOnsetMin);
minInd          = minInds(1);

StInflTime = time(minInd);
SpInflTime = time(maxInd);

adjPresOffsetMax = max(adjPresOffsetM);
[maxInds] = find(adjPresOffsetM < 0.95*adjPresOffsetMax);
maxInd    = maxInds(1);

adjPresOffsetMin = min(adjPresOffsetM);
[minInds] = find(adjPresOffsetM > 0.98*adjPresOffsetMin & adjPresOffsetM < 0.97*adjPresOffsetMin);
minInd    = minInds(1);  

fD = 20*[diff(adjPresOffsetM); 0];
mfD = mean(fD(1:200));
bumpsSt = find(fD < 50*mfD);    
bumpsSp = find(fD < 10*mfD); 

StDeflTime = time(bumpsSt(1));
SpDeflTime = time(bumpsSp(end));

adjustPres = [adjPresOnsetM, adjPresOnsetE, adjPresOffsetM, adjPresOffsetE];
InflDeflT  = [StInflTime SpInflTime StDeflTime SpDeflTime];
end

function statTable = packStatTable(ss)

varNames = {'SubjID', 'Age', 'Gender', 'f0', 'AudFB', 'StimMag', 'RespMag', 'RespPer'};
statTable = table(ss.obvSubj, ss.obvAge, ss.obvGender, ss.f0, ss.obvAudFB, ss.obvRespVar(:,2), ss.obvRespVar(:,3), ss.obvRespVar(:,4), 'VariableNames', varNames);
end

function statTable = packStatTableSingle(pooledRunStr)

numSubj = length(pooledRunStr);
ss.obvSubj = {};
ss.obvAge  = [];
ss.obvGender = {};
ss.f0        = [];
ss.obvAudFB  = {};
ss.StimMag   = [];
ss.StimMagSD = [];
ss.RespMag   = [];
ss.RespMagSD = [];
ss.RespPer   = [];
ss.RespPerSD = [];
for ii = 1:numSubj
    curRes = pooledRunStr(ii);
    numRun = length(curRes.runs);
    for jj = 1:numRun
        ss.obvSubj   = cat(1, ss.obvSubj, curRes.studyID);
        ss.obvAge    = cat(1, ss.obvAge, curRes.age);
        ss.obvGender = cat(1, ss.obvGender, curRes.gender);
        ss.f0        = cat(1, ss.f0, curRes.f0(jj));
        ss.obvAudFB  = cat(1, ss.obvAudFB, curRes.AudFB{jj}{1});
        ss.StimMag   = cat(1, ss.StimMag, curRes.respVarSingle{jj, 1}.stimMagM);
        ss.StimMagSD = cat(1, ss.StimMagSD, curRes.respVarSingle{jj, 1}.stimMagSD);
        ss.RespMag   = cat(1, ss.RespMag, curRes.respVarSingle{jj, 1}.respMagM);
        ss.RespMagSD = cat(1, ss.RespMagSD, curRes.respVarSingle{jj, 1}.respMagSD);
        ss.RespPer   = cat(1, ss.RespPer, curRes.respVarSingle{jj, 1}.respPerM);
        ss.RespPerSD = cat(1, ss.RespPerSD, curRes.respVarSingle{jj, 1}.respPerSD);   
    end
end

varNames = {'SubjID', 'Age', 'Gender', 'f0', 'AudFB', 'StimMag', 'StimMagSD', 'RespMag', 'RespMagSD', 'RespPer', 'RespPerSD',};
statTable = table(ss.obvSubj, ss.obvAge, ss.obvGender, ss.f0, ss.obvAudFB, ss.StimMag, ss.StimMagSD, ss.RespMag, ss.RespMagSD, ss.RespPer, ss.RespPerSD,'VariableNames', varNames);
end

function organizeAndPrintDemographicStats(dirs, allSubjRes)

pAnalysis = allSubjRes.pAnalysis;
ages    = allSubjRes.age;
genders = allSubjRes.gender;

meanAge  = round(mean(ages), 1);
sdAge    = round(std(ages), 1);
minAge   = round(min(ages), 1);
maxAge   = round(max(ages), 1);
rangeAge = [minAge maxAge];

numMales   = sum(strcmp(genders, 'male'));
numFemales = sum(strcmp(genders, 'female'));

genderRatio = [numMales numFemales];

maleAges   = ages(strcmp(genders, 'male'));
femaleAges = ages(strcmp(genders, 'female'));
genderAges = [maleAges, femaleAges];
gaGRP = [zeros(numMales, 1); ones(numFemales, 1)];
boxPlotLabel = {['Males (n=' num2str(numMales) ')'], ['Females (n=' num2str(numFemales) ')']};

measBox = figure('Color', [1 1 1]);
plotpos = [30 30]; plotdim = [700 1000];
set(measBox, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

fontN = 'Arial';
axisLSize = 25;

boxplot(genderAges, gaGRP, 'Labels', boxPlotLabel)
xlabel('Gender')
ylabel('Age (years)')
box off

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold')

dirs.BoxPlotFigureFile = fullfile(dirs.SavResultsDir, [pAnalysis 'AgeGenderBoxPlot.jpg']);
export_fig(dirs.BoxPlotFigureFile)

fprintf('\nSubjects in this data set are between the ages of %.1f and %.1f years (Mean: %.1f years, SD: %.1f years)\n', rangeAge(1), rangeAge(2), meanAge, sdAge)
fprintf('This data set includes %d males, and %d females\n', genderRatio(1), genderRatio(2))
end

function organizeAndSaveExcludedTrialTable(dirs, pA, allSubjRes, tossTrialTracker, tVN, svFile)

fprintf('\nAcross all subjects, %d trials were excluded Automatically:\n', allSubjRes.tossedAll);
fprintf('%d trials due to late starts\n', allSubjRes.tossedLate);
fprintf('%d trials due to voice breaks\n', allSubjRes.tossedBreak);
fprintf('%d trials due to pitch miscalc\n', allSubjRes.tossedMisCalc);
fprintf('An additional %d trials were excluded Manually, of the %d trials Manually selected\n', allSubjRes.tossedAutoMiss, allSubjRes.tossedManual);
fprintf('Automatic trial exclusion methods accounted for %s%% of trials selected Manually\n', allSubjRes.autoSuccessPerc)
fprintf('\n')
tossedTable = table(tossTrialTracker.curSess,... 
                    tossTrialTracker.tossedTrials,...
                    tossTrialTracker.manuallyExcl,...
                    tossTrialTracker.autoMiss,...
                    tossTrialTracker.perCaught,...
                    'VariableNames', tVN);
if svFile == 1                
    uitable('Data', tossedTable{:,:},...
            'ColumnName', tossedTable.Properties.VariableNames,...
            'Units', 'Normalized',...
            'Position', [0, 0, 1, 1]);
        
    dirs.excludedTrialTable = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'ExcludedTrial.csv']);
    writetable(tossedTable, dirs.excludedTrialTable, 'WriteVariableNames',true)
end
end