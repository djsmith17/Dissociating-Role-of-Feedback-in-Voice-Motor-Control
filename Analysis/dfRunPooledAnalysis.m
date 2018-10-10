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
% pooled data set at line 17.
% 
% Requires the Signal Processing Toolbox

close all
pA.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control'; 
pA.pAnalysis     = 'MaskingStudy'; % Change this name to load different pooled data sets Ex: SfN2017, LarynxPos

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
            disp('ERROR: The Results file for Run %s does not exist yet\n', run)
            return
        else   
            load(dirs.SavFile)
            % Returns a results struture of 'res'
        end
        
        subjRes = cat(2, subjRes, res);  
    end
    allDataStr = cat(1, allDataStr, subjRes);
end

allSubjRes         = initSortedStruct(pA.numCond);
allSubjRes.subject = 'Mean Participant Response';
allSubjRes.curSess = allSubjRes.subject;
allSubjRes.cond    = pA.cond;

[tossTrialTracker, tVN] = initTossedTrialTracker();

for ii = 1:pA.numPart
    participant = pA.participants{ii};
    fprintf('Sorting task conditions for %s\n', participant)
    
    sortStruc         = initSortedStruct(pA.numCond);
    sortStruc.subject = ['Participant ' num2str(ii)]; % Pooled Analysis Name
    
    sortStruc.curSess = sortStruc.subject;
    sortStruc.cond    = pA.cond;
 
    for jj = 1:pA.numRuns
        curRes = allDataStr(ii, jj);

        sortStruc.studyID = curRes.subject; % Study ID
        sortStruc.gender  = curRes.gender;
        sortStruc.age     = curRes.age;
        
        [tossCounts, tossTrialTracker, autoMiss] = combineTossedTrialTracker(curRes, tossTrialTracker);
        
        if ~isempty(autoMiss)
            curRes = adjustCurRes(curRes, autoMiss);
        end
        
        sortStruc  = combineCondTrials(pA, curRes, sortStruc, tossCounts);       
        allSubjRes = combineCondTrials(pA, curRes, allSubjRes, tossCounts);
    end
        
    sortStruc = meanCondTrials(pA, sortStruc);
    sortStruc.pltName = pA.pltNameMVi{ii};
    
    pooledRunStr(ii)   = sortStruc;
    
    allSubjRes.gender{ii} = sortStruc.gender;
    allSubjRes.age(ii)    = sortStruc.age;
end

allSubjRes = meanCondTrials(pA, allSubjRes);
allSubjRes.pltName  = pA.pltNameMVm;

fprintf('\nAcross all subjects, %d trials were excluded Automatically:\n', allSubjRes.tossedAll);
fprintf('%d trials due to late starts\n', allSubjRes.tossedLate);
fprintf('%d trials due to voice breaks\n', allSubjRes.tossedBreak);
fprintf('%d trials due to pitch miscalc\n', allSubjRes.tossedMisCalc);
fprintf('An additional %d trials were excluded Manually, of the %d trials Manually selected\n', allSubjRes.tossedAutoMiss, allSubjRes.tossedManual);
fprintf('Automatic trial exclusion methods accounted for %s%% of trials selected Manually\n', allSubjRes.autoSuccessPerc)

tossedTable = table(tossTrialTracker.curSess,... 
                    tossTrialTracker.tossedTrials,...
                    tossTrialTracker.manuallyExcl,...
                    tossTrialTracker.autoMiss,...
                    tossTrialTracker.perCaught,...
                    'VariableNames', tVN);
                
uitable('Data', tossedTable{:,:},...
        'ColumnName', tossedTable.Properties.VariableNames,...
        'Units', 'Normalized',...
        'Position', [0, 0, 1, 1]);
    
allSubjStatTable = displayStats(allSubjRes);    

[mAge, rAge, gRatio] = demoStats(allSubjRes);
fprintf('\nSubjects in this data set are between the ages of %.1f and %.1f (Mean: %.1f)\n', rAge(1), rAge(2), mAge)
fprintf('This data set includes %d males, and %d females\n\n', gRatio(1), gRatio(2))

% Save the Pooled Results
dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'ResultsDRF.mat']);
fprintf('Saving Pooled Analysis for %s\n', pA.pAnalysis)
save(dirs.SavResultsFile, 'pooledRunStr', 'allSubjRes')

dirs.behavioralResultTable = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'BehavioralResultTable.xlsx']);
writetable(allSubjStatTable, dirs.behavioralResultTable, 'WriteVariableNames',true)

dirs.excludedTrialTable = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'ExcludedTrial.xlsx']);
writetable(tossedTable, dirs.excludedTrialTable, 'WriteVariableNames',true)
end

function sortStr = initSortedStruct(numCond)
% sortStr = initSortedStruct(numCond) initializes the structure that will
% store the pooled results for each subject, or group of subjects. It is
% created to have different sizes, based on the number of conditions that
% are being tested against. I think this should generalize to subconditions
% of conditions, or two condition crossing, but I have not tested that, and
% currently the above scripts only consider one condition to test against. 

% Basic info about the session, the recordings, the subjects
sortStr.subject = [];
sortStr.gender  = [];
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

sortStr.secTimeP        = [];
sortStr.sensorPSec      = [];

sortStr.audioMf0MeanPert = cell(numCond, 1);
sortStr.audioMf0MeanCont = [];
sortStr.sensorPMean      = [];

sortStr.tossedAll        = 0;
sortStr.tossedLate       = 0;
sortStr.tossedBreak      = 0;
sortStr.tossedMisCalc    = 0;
sortStr.tossedManual     = 0;
sortStr.tossedAutoMiss   = 0;

sortStr.respVar          = cell(numCond, 1);
sortStr.respVarM         = cell(numCond, 1);

sortStr.obvSubj          = [];
sortStr.obvAge           = [];
sortStr.obvGender        = [];
sortStr.obvAudFB         = [];
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

whichCondAr = strcmp(pA.cond, eval(pA.condVar));
wC          = find(whichCondAr == 1);            % Which Condition?

polRes.runs{wC}   = cat(1, polRes.runs{wC}, {curRes.run});
polRes.runf0b{wC} = cat(1, polRes.runf0b{wC}, curRes.f0b);

polRes.AudFB{wC}  = cat(1, polRes.AudFB{wC}, {curRes.AudFB});

polRes.allContTrials     = cat(1, polRes.allContTrials, curRes.numContTrialsFin);
polRes.allPertTrials{wC} = cat(1, polRes.allPertTrials{wC}, curRes.numPertTrialsFin);

polRes.secTime             = curRes.secTime;
polRes.audioMf0SecPert{wC} = cat(2, polRes.audioMf0SecPert{wC}, curRes.audioMf0SecPert);
polRes.audioMf0SecCont     = cat(2, polRes.audioMf0SecCont, curRes.audioMf0SecCont);
polRes.respVar{wC}         = cat(1, polRes.respVar{wC}, curRes.respVar);

polRes.secTimeP            = curRes.secTimeP;
polRes.sensorPSec          = cat(2, polRes.sensorPSec, curRes.sensorPSec);

polRes.tossedAll      = polRes.tossedAll + tossT.A;       % Total Automatic Excluded Trials
polRes.tossedLate     = polRes.tossedLate + tossT.L;      % Late Start
polRes.tossedBreak    = polRes.tossedBreak + tossT.B;     % Voice Break
polRes.tossedMisCalc  = polRes.tossedMisCalc + tossT.C;   % f0 Miscalc
polRes.tossedManual   = polRes.tossedManual + tossT.M;    % Total Manual Excluded Trials
polRes.tossedAutoMiss = polRes.tossedAutoMiss + tossT.aM; % Trials Manually removed, but missed by auto methods.

[respVar, ~, ~]      = InflationResponse(curRes.secTime, curRes.audioMf0SecPert);
polRes.respVar{wC}   = respVar;

ages = repmat(curRes.age, curRes.numPertTrialsFin, 1);
subjs = cell(1, curRes.numPertTrialsFin);
[subjs{:}] = deal(curRes.subject);
genders = cell(1, curRes.numPertTrialsFin);
[genders{:}] = deal(curRes.gender);
AudFBs = cell(1, curRes.numPertTrialsFin);
[AudFBs{:}] = deal(curRes.AudFB);

polRes.obvSubj         = cat(1, polRes.obvSubj, subjs');
polRes.obvAge          = cat(1, polRes.obvAge, ages);
polRes.obvGender       = cat(1, polRes.obvGender, genders');
polRes.obvAudFB        = cat(1, polRes.obvAudFB, AudFBs');
polRes.obvRespVar      = cat(1, polRes.obvRespVar, respVar);
end

function [tossCounts, tossTrialTracker, autoMiss] = combineTossedTrialTracker(curRes, tossTrialTracker)

[tossCounts, aTossTrials, aTossDetails] = identifyTossedTrials(curRes);

[mTossTrials, mTossDetails] = unpackManuallyExcludedTrials(curRes.incTrialInfo);

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

acceptedContTrials = curRes.contIdxFin;
[~, contKeepIdx] = setdiff(acceptedContTrials, autoMiss);
numCont2Keep = length(contKeepIdx);

acceptedPertTrials = curRes.pertIdxFin;
[~, pertKeepIdx] = setdiff(acceptedPertTrials, autoMiss);
numPert2Keep = length(pertKeepIdx);

curRes.numContTrialsFin = numCont2Keep;
curRes.numPertTrialsFin = numPert2Keep;

curRes.audioMf0SecCont = curRes.audioMf0SecCont(:, contKeepIdx, :);
curRes.audioMf0SecPert = curRes.audioMf0SecPert(:, pertKeepIdx, :);
curRes.sensorPSec      = curRes.sensorPSec(:, pertKeepIdx, :);
end

function polRes = meanCondTrials(pA, polRes)

polRes.numContTrialsFin = sum(polRes.allContTrials);
polRes.audioMf0MeanCont = meanSecData(polRes.audioMf0SecCont);
polRes.sensorPMean      = meanSecData(polRes.sensorPSec);

for wC = 1:pA.numCond
    polRes.f0b(wC)              = mean(polRes.runf0b{wC});
    
    polRes.numPertTrialsFin(wC) = sum(polRes.allPertTrials{wC});
    polRes.audioMf0MeanPert{wC} = meanSecData(polRes.audioMf0SecPert{wC});
    polRes.respVarM{wC}         = mean(polRes.respVar{wC}, 1);
end

% Identify Limits of the newly meaned data
lims = identifyLimits(polRes);
polRes.limitsAmean = lims.audioMean;
polRes.limitsPmean = lims.presMean;

% Scale Pressure traces against the f0 traces
[sensorPAdjust, InflDeflT] = adjustPressureVals(polRes);
polRes.sensorPAdjust = sensorPAdjust;
polRes.InflDeflT     = InflDeflT;

% Identify success of automatic trial exclusion methods
if polRes.tossedManual > 0
    perc = round(100*(1-(polRes.tossedAutoMiss/polRes.tossedManual)));
    polRes.autoSuccessPerc = num2str(perc);
else
    polRes.autoSuccessPerc = 'N/A';
end

statTable        = packStatTable(polRes);
polRes.statTable = statTable;
end

function meanData = meanSecData(secData)

OnsetSecs  = secData(:,:,1);
OffsetSecs = secData(:,:,2);
[~, numTrial] = size(OnsetSecs);

meanOnset  = mean(OnsetSecs, 2);
meanOffset = mean(OffsetSecs, 2);

stdOnset   = std(OnsetSecs, 0, 2);
stdOffset  = std(OffsetSecs, 0, 2);

SEMOnset   = stdOnset/sqrt(numTrial);  % Standard Error
SEMOffset  = stdOffset/sqrt(numTrial); % Standard Error

NCIOnset   = 1.96*SEMOnset;  % 95% Confidence Interval
NCIOffset  = 1.96*SEMOffset; % 95% Confidence Interval

meanData   = [meanOnset NCIOnset meanOffset NCIOffset];
end

function [respVar, respVarM, respVarSD] = InflationResponse(secTime, secAudio)
% [respVar, respVarm, respVarSD, InflaStimVar] = InflationResponse(secTime, secAudio)
% Identifies the relevant pitch contour characteristics that are important
% for deciding how a participant responded to the inflation of the balloon
% during production. iR is a structure representing the result variables
% from studying the inflation response (iR). 
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

[numSamp, numTrial, ~] = size(secAudio); % Size of the data we are dealing with

ir.numSamp  = numSamp;
ir.numTrial = numTrial;
ir.time     = secTime;

ir.iAtOnset = find(secTime == 0); % Index where t = 0
ir.tAtOnset = 0;                  % Time at t = 0 ...duh
ir.vAtOnset = [];                 % f0 value at t = 0

ir.iPostOnsetR = find(0 <= secTime & .20 >= secTime); % Range of indices between t = 0ms and t = 200ms;
ir.iAtMin  = [];                  % Index at min f0 value in PostOnsetR
ir.tAtMin  = [];                  % Time at min f0 value in PostOnsetR
ir.vAtMin  = [];                  % Min f0 value in PostOnsetR
ir.stimMag = [];                  % ir.vAtMin - ir.vAtOnset ..in a perfect world vAtOnset = 0

ir.iAtResp = ir.numSamp;          % Index of f0 value when participant 'fully' responded...right now = last value in section
ir.tAtResp = ir.time(ir.numSamp); % Time at f0 value when participant 'fully' responded
ir.vAtResp = [];                  % f0 value when participant 'fully' responded 
ir.respMag = [];                  % vAtResp - vAtMin   ...distance traveled
ir.respPer = [];                  % Percent change from stimMag to respMag

% Variables to be concatenated and saved as outputs 
tAtMin  = []; stimMag = [];
respMag = []; respPer = [];
for i = 1:numTrial
    onset = secAudio(:,i,1); % Go trial by trial; First 3D layer is Onset
    ir.vAtOnset = onset(ir.iAtOnset); % f0 value at t = 0

    [minOn, minIdx] = min(onset(ir.iPostOnsetR)); % Minimum f0 in PostOnsetR
    ir.iAtMin = ir.iPostOnsetR(minIdx);           % Indice of the min f0 value
    ir.tAtMin = ir.time(ir.iAtMin);               % Time at min f0 value in PostOnsetR
    ir.vAtMin = minOn;                            % Min f0 value in PostOnsetR
    ir.stimMag = ir.vAtMin - ir.vAtOnset;         % Distance traveled from onset to min value
    
    ir.vAtResp = onset(ir.iAtResp);               % f0 value when participant 'fully' responded 
    ir.respMag = ir.vAtResp - ir.vAtMin;          % Distance traveled from min f0 value to response f0 value
    
    ir.respPer = 100*(ir.respMag/abs(ir.stimMag));% Percent change from stimMag to respMag 
    
    if ir.stimMag == 0
        ir.respPer = 0.0;
    end
    
    % Concatenate the results from this trial 
    tAtMin   = cat(1, tAtMin, ir.tAtMin);
    stimMag  = cat(1, stimMag, ir.stimMag); 
    respMag  = cat(1, respMag, ir.respMag); 
    respPer  = cat(1, respPer, ir.respPer);
end

% Organize the results 
respVar   = [tAtMin stimMag respMag respPer];
respVarM  = mean(respVar, 1);
respVarSD = std(respVar, 0, 1);
end

function lims = identifyLimits(ss)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f0 Traces Limits
mf0MeanPert = ss.audioMf0MeanPert;
numCond     = length(mf0MeanPert);

setupBoundSec = zeros(numCond, 1);
setlwBoundSec = zeros(numCond, 1);
for ii = 1:numCond
    audioMean = mf0MeanPert{ii};    

    [~, Imax] = max(audioMean(:,1)); %Max Pert Onset
    upBoundOn = round(audioMean(Imax,1) + audioMean(Imax,2) + 10);
    [~, Imin] = min(audioMean(:,1)); %Min Pert Onset
    lwBoundOn = round(audioMean(Imin,1) - audioMean(Imin,2) - 10);

    [~, Imax] = max(audioMean(:,3)); %Max Pert Offset
    upBoundOf = round(audioMean(Imax,3) + audioMean(Imax,4) + 10);
    [~, Imin] = min(audioMean(:,3)); %Min Pert Offset
    lwBoundOf = round(audioMean(Imin,3) - audioMean(Imin,4) - 10);

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
    
    setupBoundSec(ii) = upBoundSec;
    setlwBoundSec(ii) = lwBoundSec;  
end

maxUpBound = max(setupBoundSec); % Max f0 Bound
minLwBound = min(setlwBoundSec); % Min f0 Bound 

lims.audioMean = [-0.5 1.0 minLwBound maxUpBound];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pressure Traces Limits
maxPres = max(ss.sensorPMean(:,1)) + 0.5;
minPres = min(ss.sensorPMean(:,1)) - 0.1;

lims.presMean = [-0.5 1.0 minPres maxPres];
end

function [adjustPres, InflDeflT] = adjustPressureVals(ss)

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

function statLib = packStatLib(ss)

cond1 = ss.respVar{1};
cond2 = ss.respVar{2};

condM1 = ss.respVarM(1,:);
condM2 = ss.respVarM(2,:);

[~, pStim] = ttest2(cond1(:,2), cond2(:,2));
[~, pResp] = ttest2(cond1(:,3), cond2(:,3));
[~, pPerc] = ttest2(cond1(:,4), cond2(:,4));

statLib(1) = condM1(2); % Condition 1 StimMag
statLib(2) = condM2(2); % Condition 2 StimMag
statLib(3) = condM1(3); % Condition 1 RespMag
statLib(4) = condM2(3); % Condition 2 RespMag
statLib(5) = condM1(4); % Condition 1 %
statLib(6) = condM2(4); % Condition 2 %
statLib(7) = pStim;     % p-value stimulus
statLib(8) = pResp;     % p-value response
statLib(9) = pPerc;     % p-value percent increase 
end

function statTable = packStatTable(ss)

varNames = {'SubjID', 'Age', 'Gender', 'AudFB', 'StimMag', 'RespMag', 'RespPer'};
statTable = table(ss.obvSubj, ss.obvAge, ss.obvGender, ss.obvAudFB, ss.obvRespVar(:,2), ss.obvRespVar(:,3), ss.obvRespVar(:,4), 'VariableNames', varNames);
end

function [meanAge, rangeAge, genderRatio] = demoStats(allSubjRes)

ages    = allSubjRes.age;
genders = allSubjRes.gender;

meanAge = round(mean(ages), 1);
minAge  = round(min(ages), 1);
maxAge  = round(max(ages), 1);
rangeAge = [minAge maxAge];

numMales = sum(strcmp(genders, 'male'));
numFemales = sum(strcmp(genders, 'female'));

genderRatio = [numMales numFemales];
end

function allSubjStatTable = displayStats(allSubjRes)

allSubjStatTable = allSubjRes.statTable;
allSM = allSubjStatTable.StimMag;
allRM = allSubjStatTable.RespMag;
allRP = allSubjStatTable.RespPer;
aSRM = fitrm(allSubjStatTable, 'StimMag-RespPer~AudFB');
        
normDist = figure('Color', [1 1 1]);
subplot(1,3,1)
histogram(allSM)
title('Stimulus Magnitude')
box off

subplot(1,3,2)
histogram(allRM)
title('Response Magnitude')
box off

subplot(1,3,3)
histogram(allRP)
title('Response Percentage')
box off



end