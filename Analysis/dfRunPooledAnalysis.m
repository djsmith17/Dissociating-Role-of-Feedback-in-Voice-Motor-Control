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

profile on

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

% Extract variables from the configuration file
pA.participants  = cF.participants; % List of multiple participants.
pA.numPart       = length(pA.participants);
pA.runs          = cF.runs;         % All runs to consider 
pA.numRuns       = length(pA.runs);
pA.totnumRuns    = pA.numPart*pA.numRuns; % How many total runs do we care about?
pA.cond          = cF.cond;         % Conditions to test against
pA.numCond       = length(pA.cond); 
pA.condVar       = cF.condVar;      % Variable to test the condition
pA.testExt       = cF.testExt;
pA.pubCond       = cF.pubCond;
clear cF

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
clear cF participant run subjRes res

numObs             = pA.numPart*pA.numCond;
allSubjRes         = initSortedStruct(pA, numObs);
allSubjRes.subject = 'Mean Participant Response';
allSubjRes.curSess = allSubjRes.subject;
allSubjRes.cond    = pA.cond;
allSubjRes.pubCond = pA.pubCond;
allSubjRes.statTableSingle = initStatTableSingle(numObs);

tossTrialTable = initTossedTrialTable(pA.totnumRuns);

%Sort all of the loaded analyzed data
runItr = 0;
for ii = 1:pA.numPart
    participant = pA.participants{ii};
    fprintf('Sorting task conditions for %s\n', participant)
    
    sortStr         = initSortedStruct(pA, 2);
    sortStr.subject = ['Participant ' num2str(ii)]; % Pooled Analysis Name
    
    sortStr.curSess = sortStr.subject;
    sortStr.cond    = pA.cond;
    sortStr.pubCond = pA.pubCond;
 
    wcCount_indivi = [0 0];
    wcCount_group  = [0 0];
    for jj = 1:pA.numRuns
        curRes = allDataStr(ii, jj);

        sortStr.studyID = curRes.subject; % Study ID
        sortStr.expType = curRes.expType;
        sortStr.gender  = curRes.gender;
        sortStr.age     = round(curRes.age, 1);
        
        runItr = runItr + 1;
        [tossCounts, tossTrialTable, autoMiss] = combineTossedTrialTracker(curRes, tossTrialTable, runItr);
        
        if ~isempty(autoMiss)
            curRes = adjustCurRes(curRes, autoMiss);
        end
        
        [sortStr, wcCount_indivi] = combineCondTrials(pA, curRes, sortStr, tossCounts, wcCount_indivi);       
        [allSubjRes, wcCount_group] = combineCondTrials(pA, curRes, allSubjRes, tossCounts, wcCount_group);
    end
        
    sortStr = meanCondTrials(pA, sortStr);
    sortStr.pltName = pA.pltNameMVi{ii};
    
    % Organize Stat Tables
    sortStr  = popStatTableObs(pA, sortStr, sortStr, 1);
    allSubjRes = popStatTableObs(pA, allSubjRes, sortStr, ii);
    allSubjRes = packStatTableSingle(pA, allSubjRes, sortStr, ii);
    
    pooledRunStr(ii)   = sortStr;
    
    allSubjRes.expType    = sortStr.expType;
    allSubjRes.gender{ii} = sortStr.gender;
    allSubjRes.age(ii)    = sortStr.age;
end
clear ii jj participant runItr numObs curRes

allSubjRes = catSubjMeans(pA, allSubjRes, pooledRunStr);
allSubjRes = meanCondTrials(pA, allSubjRes);
allSubjRes.pltName = pA.pltNameMVm;

% Organize and Print the Stats of the Demographics included in this study
organizeAndPrintDemographicStats(dirs, allSubjRes);

% Organize and Save the Table of Excluded Trials
organizeAndSaveExcludedTrialTable(dirs, pA, allSubjRes, tossTrialTable, 1)

organizePooledResultsForFrank(dirs, allSubjRes)

% Save the Pooled Results
dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'ResultsDRF.mat']);
fprintf('Saving Pooled Analysis for %s\n', pA.pAnalysis)
save(dirs.SavResultsFile, 'pooledRunStr', 'allSubjRes')

%Apply Stats as appropriate
if strcmp(pA.pAnalysis, 'MaskingStudy')
    drawExpPressureDist(dirs, pA, pooledRunStr)
    StatsOrg_MaskingNoiseStudy(dirs, pA, allSubjRes);
elseif strcmp(pA.pAnalysis, 'DRF_Som')
    drawSubjRespVarVoicingCorr(dirs, pA, pooledRunStr)
    drawExpPressureDist(dirs, pA, pooledRunStr)
    StatsOrg_DRF_Som(dirs, pA, allSubjRes);
elseif strcmp(pA.pAnalysis, 'DRF_Aud')
    drawSubjRespVarVoicingCorr(dirs, pA, pooledRunStr)
    drawSubjRespVarDists_AllTrial(dirs, pooledRunStr)
    StatsOrg_DRF_Aud(dirs, pA, allSubjRes);
end

profile viewer
profsave(profile('info'), fullfile(dirs.SavResultsDir, 'Improved Remove Trial'))
end

function sortStr = initSortedStruct(pA, numObs)
% sortStr = initSortedStruct(pA, numObs) initializes the structure that will
% store the sorted, pooled results for each subject, or group of subjects. It is
% created to have different sizes, based on the number of conditions that
% are being tested against. I think this should generalize to subconditions
% of conditions, or two condition crossing, but I have not tested that, and
% currently the above scripts only consider one condition to test against. 

pAnalysis = pA.pAnalysis;
numCond = pA.numCond;

% Basic info about the session, the recordings, the subjects
sortStr.pAnalysis = pAnalysis; % Pooled Analysis Name.
sortStr.expType = [];          % Type of experiment for these data
sortStr.studyID = [];          % Study ID of the participant
sortStr.subject = [];          % Numbered participant name. Separate from StudyID
sortStr.gender  = [];          % participant gender
sortStr.age     = [];          % participant age
sortStr.curSess = [];          % current session being analyzed. Currently replicates subject field and nothing more.

sortStr.runs    = cell(numCond, 1);  % 
sortStr.cond    = cell(numCond, 1);
sortStr.pubCond = cell(numCond, 1);
sortStr.runf0b  = cell(numCond, 1);
sortStr.f0b     = zeros(numCond, 1);

sortStr.AudFB   = cell(numCond, 1);

% Trial Information
sortStr.allContTrials    = [];
sortStr.numContTrialsFin = 0;
sortStr.allPertTrials    = cell(numCond, 1);
sortStr.numPertTrialsFin = zeros(numCond, 1);
sortStr.AuMHDelays       = cell(numCond, 1);
sortStr.AuMHDelayMean    = [];

% Audio Dynamics
sortStr.secTime         = [];
sortStr.audioMf0SecPert = cell(numCond, 1);
sortStr.audioMf0SecCont = [];
sortStr.audioHf0SecPert = cell(numCond, 1);
sortStr.audioHf0SecCont = [];

% Pre Voicing Dynamics
sortStr.prePertVoicingTimeCont = [];
sortStr.prePertVoicingTimePert = cell(numCond,1);

% Sensor Dynamics
sortStr.secTimeP        = [];
sortStr.sensorPSec      = [];
sortStr.sensorPOnOff    = [];
sortStr.sensorPRiseTime = [];

sortStr.audioMf0MeanPert = cell(numCond, 1);
sortStr.audioMf0MeanCont = [];
sortStr.audioHf0MeanPert = cell(numCond, 1);
sortStr.audioHf0MeanCont = [];
sortStr.sensorPMean      = [];
sortStr.sensorPOnOffMean = [];
sortStr.sensorPRiseTimeM = [];

% Trial Tossing Record
sortStr.tossedAll        = 0;
sortStr.tossedLate       = 0;
sortStr.tossedBreak      = 0;
sortStr.tossedMisCalc    = 0;
sortStr.tossedManual     = 0;
sortStr.tossedAutoMiss   = 0;

sortStr.pertLengths    = cell(numCond, 1);
sortStr.pertOnsetTimes = cell(numCond, 1);

sortStr.respVarSingle    = cell(numCond, 1);
sortStr.respVarM         = cell(numCond, 1);

sortStr.statTable = initStatTable(numObs);

sortStr.audioMf0SecPert_all  = cell(numCond, 2);
sortStr.audioHf0SecPert_all  = cell(numCond, 2);
sortStr.audioMf0MeanPert_all = cell(numCond, 2);
sortStr.audioHf0MeanPert_all = cell(numCond, 2);
sortStr.respVarM_all         = cell(numCond, 2);
end

function [tossTrialTable] = initTossedTrialTable(totnumRuns)
% Initalizes the table which organizes which trials were tossed, how they
% were tossed, and why they were tossed.

genVar = cell(totnumRuns, 1);

% Tossed Trials Table Variable Names
tVN = {'CurSess', 'AutoRemoved', 'ManuallyRemoved', 'MissedByAuto', 'PercentCaught', 'TotalRemoved', 'TotalTrials', 'PercentRemoved'};

tossTrialTable = table(genVar, genVar, genVar, genVar, genVar, genVar, genVar, genVar);
tossTrialTable.Properties.VariableNames = tVN;
end

function [statTable] = initStatTable(numObs)
%Initalize the stat table which details the basic outcomes of the sorting
%and pooling of analyzed results. Specifically this organizes the outcome
%variables of tAtMin, StimMag, RespMag, RespPer and f0 measured from the 
%mean  run f0-traces from each participant. 

varNames = {'SubjID', 'Age', 'Gender', 'f0', 'AudFB', 'tAtMin', 'StimMag', 'RespMag', 'RespPer', 'tAtMin1', 'tAtMin2', 'StimMag1', 'StimMag2', 'RespMag1', 'RespMag2', 'RespPer1', 'RespPer2'};
varTypes = {'string' 'double', 'string', 'double', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
numVar = length(varNames);

statTable = table('Size', [numObs numVar], 'VariableTypes', varTypes, 'VariableNames', varNames);
end

function [statTableSingle] = initStatTableSingle(numObs)
%Initalize the stat table which details the basic outcomes of the sorting
%and pooling of analyzed results. Specifically this organizes the outcome
%variables of tAtMin, StimMag, RespMag, RespPer and f0 measured from the 
%mean  run f0-traces from each participant. 

varNames = {'SubjID', 'Age', 'Gender', 'f0', 'AudFB', 'StimMag', 'StimMagSD', 'RespMag', 'RespMagSD', 'RespPer', 'RespPerSD'};
varTypes = {'string' 'double', 'string', 'double', 'string', 'double', 'double', 'double', 'double', 'double', 'double'};
numVar = length(varNames);

statTableSingle = table('Size', [numObs numVar], 'VariableTypes', varTypes, 'VariableNames', varNames);
end

function [sortStr, wCCount] = combineCondTrials(pA, curRes, sortStr, tossT, wCCount)
%combineCondTrials organizes the information that were sorted from each run
%and combines them into  

AD = curRes.audioDynamics; % Audio Dynamics
PD = curRes.presSDsv;      % Pressure Dynamics

whichCondAr = strcmp(pA.cond, eval(pA.condVar));
wC          = find(whichCondAr == 1);            % Which Condition?
wCCount     = wCCount + whichCondAr; % How many times have we had this condition? 

% Run Information
sortStr.runs{wC}   = cat(1, sortStr.runs{wC}, {curRes.run});
sortStr.runf0b{wC} = cat(1, sortStr.runf0b{wC}, curRes.f0b);

sortStr.AudFB{wC}  = cat(1, sortStr.AudFB{wC}, {curRes.AudFB});

% Trial Information
sortStr.allContTrials     = cat(1, sortStr.allContTrials, curRes.numContTrialsFin);
sortStr.allPertTrials{wC} = cat(1, sortStr.allPertTrials{wC}, curRes.numPertTrialsFin);
sortStr.AuMHDelays{wC}    = cat(1, sortStr.AuMHDelays{wC}, curRes.AuMHDelaysinc);

% Audio Dynamics
sortStr.secTime             = curRes.secTime;
sortStr.audioMf0SecPert{wC} = cat(2, sortStr.audioMf0SecPert{wC}, curRes.audioMf0SecPert);
sortStr.audioHf0SecPert{wC} = cat(2, sortStr.audioHf0SecPert{wC}, curRes.audioHf0SecPert);
sortStr.audioMf0SecCont     = cat(2, sortStr.audioMf0SecCont, curRes.audioMf0SecCont);
sortStr.audioHf0SecCont     = cat(2, sortStr.audioHf0SecCont, curRes.audioHf0SecCont);

% the all...harumph
sortStr.audioMf0SecPert_all{wC, wCCount(wC)} = cat(2, sortStr.audioMf0SecPert_all{wC, wCCount(wC)}, curRes.audioMf0SecPert);
sortStr.audioHf0SecPert_all{wC, wCCount(wC)} = cat(2, sortStr.audioHf0SecPert_all{wC, wCCount(wC)}, curRes.audioHf0SecPert);

% Pre Voicing Dynamics
sortStr.prePertVoicingTimeCont     = cat(1, sortStr.prePertVoicingTimeCont, curRes.prePertVoicingTimeinc(curRes.contIdxFin));
sortStr.prePertVoicingTimePert{wC} = cat(1, sortStr.prePertVoicingTimePert{wC}, curRes.prePertVoicingTimeinc(curRes.pertIdxFin));

% Sensor Dynamics
sortStr.secTimeP            = PD.timeSec;
sortStr.sensorPSec          = cat(2, sortStr.sensorPSec, PD.sensorSec);
sortStr.sensorPOnOff        = cat(1, sortStr.sensorPOnOff, PD.OnOffVal);
sortStr.sensorPRiseTime     = cat(1, sortStr.sensorPRiseTime, PD.riseTimeM);

% Trial Tossing Record
sortStr.tossedAll      = sortStr.tossedAll + tossT.A;       % Total Automatic Excluded Trials
sortStr.tossedLate     = sortStr.tossedLate + tossT.L;      % Late Start
sortStr.tossedBreak    = sortStr.tossedBreak + tossT.B;     % Voice Break
sortStr.tossedMisCalc  = sortStr.tossedMisCalc + tossT.C;   % f0 Miscalc
sortStr.tossedManual   = sortStr.tossedManual + tossT.M;    % Total Manual Excluded Trials
sortStr.tossedAutoMiss = sortStr.tossedAutoMiss + tossT.aM; % Trials Manually removed, but missed by auto methods.

% Analysis of Perturbation Lengths
if strcmp(sortStr.expType, 'Somatosensory Perturbation_Perceptual')
    pertLengths = PD.pertTime(:,2) - PD.pertTime(:,1);
    sortStr.pertLengths{wC} = cat(1, sortStr.pertLengths{wC}, pertLengths);
    
    pertOnsetTime = 0.8 + curRes.presSD.pertTime(:,1);
    sortStr.pertOnsetTimes{wC} = cat(1, sortStr.pertOnsetTimes{wC}, pertOnsetTime);
else
    sortStr.pertLengths{wC}    = 0;
    sortStr.pertOnsetTimes{wC} = 0;
end
end

function [tossCounts, tossTrialTable, autoMiss] = combineTossedTrialTracker(curRes, tossTrialTable, runItr)
% Combine information about which trials were tossed

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

%Correct the Total
tossCounts.A = tossCounts.A + tossCounts.aM;

%Percent of Trials Removed from this run
perRem = round(100*((tossCounts.A)/(curRes.numTrial)),1);

tossTrialTable.CurSess(runItr)          = {curRes.curSess};
tossTrialTable.AutoRemoved(runItr)      = {aTossDetails};
tossTrialTable.ManuallyRemoved(runItr)  = {mTossDetails};
tossTrialTable.MissedByAuto(runItr)     = {num2str(autoMiss)};
tossTrialTable.PercentCaught(runItr)    = {perCaught};

tossTrialTable.TotalRemoved(runItr)     = {num2str(tossCounts.A)};
tossTrialTable.TotalTrials(runItr)      = {num2str(curRes.numTrial)};
tossTrialTable.PercentRemoved(runItr)   = {num2str(perRem)};
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
%Mean the pooled results and create new cross run results, and cross
%participant results

polRes.numContTrialsFin = sum(polRes.allContTrials);
polRes.audioMf0MeanCont = meanSecData(polRes.audioMf0SecCont);
polRes.audioHf0MeanCont = meanSecData(polRes.audioHf0SecCont);

if strcmp(polRes.expType, 'Somatosensory Perturbation_Perceptual')
    polRes.sensorPMean = meanSecData(polRes.sensorPSec);
    polRes.sensorPOnOffMean = mean(polRes.sensorPOnOff, 1);
    polRes.sensorPRiseTimeM = mean(polRes.sensorPRiseTime, 1);
else
    polRes.sensorPMean.ON.mean = zeros(length(polRes.audioMf0SecPert{1}), 1);
    polRes.sensorPMean.OF.mean = zeros(length(polRes.audioMf0SecPert{1}), 1);
    polRes.sensorPOnOffMean = [0 0];
    polRes.sensorPRiseTimeM = [0 0];
end

[~, numTypes] = size(polRes.audioMf0MeanPert_all);

for wC = 1:pA.numCond
    polRes.f0b(wC)              = mean(polRes.runf0b{wC});
    
    polRes.numPertTrialsFin(wC) = sum(polRes.allPertTrials{wC});
    polRes.audioMf0MeanPert{wC} = meanSecData(polRes.audioMf0SecPert{wC});
    polRes.audioHf0MeanPert{wC} = meanSecData(polRes.audioHf0SecPert{wC});
    polRes.AuMHDelayMean(1)     = mean(polRes.AuMHDelays{wC});
    polRes.AuMHDelayMean(2)     = std(polRes.AuMHDelays{wC});
    
    for ii = 1:numTypes
        polRes.audioMf0MeanPert_all{wC, ii} = meanSecData(polRes.audioMf0SecPert_all{wC, ii});
        polRes.audioHf0MeanPert_all{wC, ii} = meanSecData(polRes.audioHf0SecPert_all{wC, ii});
    end
    
    if strcmp(polRes.expType, 'Somatosensory Perturbation_Perceptual')
        audioDynamicsAllTraces = InflationResponseSingle(polRes.secTime, polRes.audioMf0SecPert{wC});
        audioDynamicsMeanTrace = InflationResponse(polRes.secTime, polRes.audioMf0MeanPert{wC});
        
        for ii = 1:numTypes
            audioDynamicsMeanTrace_AllRun = InflationResponse(polRes.secTime, polRes.audioMf0MeanPert_all{wC, ii});
            polRes.respVarM_all{wC, ii} = audioDynamicsMeanTrace_AllRun.respVarM;
        end
    elseif strcmp(polRes.expType, 'Auditory Perturbation_Perceptual')
        audioDynamicsAllTraces = PitchShiftReflexResponseSingle(polRes.secTime, polRes.audioMf0SecPert{wC});
        audioDynamicsMeanTrace = InflationResponse(polRes.secTime, polRes.audioHf0MeanPert{wC});
        
        for ii = 1:numTypes
            audioDynamicsMeanTrace_AllRun = InflationResponse(polRes.secTime, polRes.audioHf0MeanPert_all{wC, ii});
            polRes.respVarM_all{wC, ii} = audioDynamicsMeanTrace_AllRun.respVarM;
        end
    end
    polRes.respVarSingle{wC} = audioDynamicsAllTraces;
    polRes.respVarM{wC}      = audioDynamicsMeanTrace.respVarM;
end

% Identify Limits of the newly meaned data
lims = identifyLimits(polRes);
polRes.limitsAmean  = lims.audioMean;
polRes.limitsMHmean = lims.audioMHMean;
polRes.limitsPmean  = lims.presMean;

% Scale Pressure traces against the f0 traces
if strcmp(polRes.expType, 'Somatosensory Perturbation_Perceptual')
    [sensorPAdjust, InflDeflT] = scalePressureVals(polRes);
elseif strcmp(polRes.expType, 'Auditory Perturbation_Perceptual')
    polRes.secTimeP = linspace(-0.5, 1.0, 12000);
    sig = -100*ones(size(polRes.secTimeP)); 

    timeLZero = polRes.secTimeP <= 0;
    sig(timeLZero) = 0;
    
    timeRampD = find((polRes.secTimeP >= 0 ) & (polRes.secTimeP <= 0.11));
    rampDValues = linspace(0, -100, length(timeRampD));
    sig(timeRampD) = rampDValues;
    sig = sig';
    polRes.sensorPMean.ON.mean = sig;
    polRes.sensorPMean.OF.mean = -1*(sig+100);
    
    sensorPAdjust = polRes.sensorPMean;
    InflDeflT     = [0.04 0.24 1.04 1.24];
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

function polRes = popStatTableObs(pA, polRes, sortStruc, curSubj)

for wC = 1:pA.numCond
    
    curIdx = (curSubj*pA.numCond) - pA.numCond + wC;
    polRes.statTable.SubjID(curIdx)  = sortStruc.studyID;
    polRes.statTable.Age(curIdx)     = sortStruc.age;
    polRes.statTable.Gender(curIdx)  = sortStruc.gender;
    polRes.statTable.f0(curIdx)      = sortStruc.f0b(wC);
    polRes.statTable.AudFB(curIdx)   = sortStruc.AudFB{wC}(1);
    polRes.statTable.tAtMin(curIdx)  = sortStruc.respVarM{wC}(1);
    polRes.statTable.StimMag(curIdx) = sortStruc.respVarM{wC}(2);
    polRes.statTable.RespMag(curIdx) = sortStruc.respVarM{wC}(3);
    polRes.statTable.RespPer(curIdx) = sortStruc.respVarM{wC}(4);
    
    polRes.statTable.tAtMin1(curIdx) = sortStruc.respVarM_all{wC,1}(1);
    polRes.statTable.tAtMin2(curIdx) = sortStruc.respVarM_all{wC,2}(1);
    polRes.statTable.StimMag1(curIdx) = sortStruc.respVarM_all{wC,1}(2);
    polRes.statTable.StimMag2(curIdx) = sortStruc.respVarM_all{wC,2}(2);
    polRes.statTable.RespMag1(curIdx) = sortStruc.respVarM_all{wC,1}(3);
    polRes.statTable.RespMag2(curIdx) = sortStruc.respVarM_all{wC,2}(3);
    polRes.statTable.RespPer1(curIdx) = sortStruc.respVarM_all{wC,1}(4);
    polRes.statTable.RespPer2(curIdx) = sortStruc.respVarM_all{wC,2}(4);
end
end

function allSubjRes = catSubjMeans(pA, allSubjRes, pooledRunStr)

numPartici    = length(pooledRunStr);
[~, numTypes] = size(allSubjRes.audioMf0SecPert_all);

audioMf0SecContAllSubj = [];
audioMf0SecPertAllSubj = cell(pA.numCond,1);
audioMf0SecPertAllSubj_all = cell(pA.numCond,2);
audioHf0SecContAllSubj = [];
audioHf0SecPertAllSubj = cell(pA.numCond,1);
audioHf0SecPertAllSubj_all = cell(pA.numCond,2);
count = 0;
for nP = 1:numPartici
    curSubj = pooledRunStr(nP);
    
    meanOnsetCont   = curSubj.audioMf0MeanCont.ON.mean;
    meanOffsetCont  = curSubj.audioMf0MeanCont.OF.mean;
    secCont(:,:,1) = meanOnsetCont;
    secCont(:,:,2) = meanOffsetCont;
    audioMf0SecContAllSubj = cat(2, audioMf0SecContAllSubj, secCont);
    
    meanOnsetCont   = curSubj.audioHf0MeanCont.ON.mean;
    meanOffsetCont  = curSubj.audioHf0MeanCont.OF.mean;
    secCont(:,:,1) = meanOnsetCont;
    secCont(:,:,2) = meanOffsetCont;
    audioHf0SecContAllSubj = cat(2, audioHf0SecContAllSubj, secCont);
    
    for wC = 1:pA.numCond
        meanOnsetPert   = curSubj.audioMf0MeanPert{wC}.ON.mean;
        meanOffsetPert  = curSubj.audioMf0MeanPert{wC}.OF.mean;
        secPert(:,:,1) = meanOnsetPert;
        secPert(:,:,2) = meanOffsetPert;
        audioMf0SecPertAllSubj{wC} = cat(2, audioMf0SecPertAllSubj{wC}, secPert);
        
        meanOnsetPert   = curSubj.audioHf0MeanPert{wC}.ON.mean;
        meanOffsetPert  = curSubj.audioHf0MeanPert{wC}.OF.mean;
        secPert(:,:,1) = meanOnsetPert;
        secPert(:,:,2) = meanOffsetPert;
        audioHf0SecPertAllSubj{wC} = cat(2, audioHf0SecPertAllSubj{wC}, secPert);
        
        for ii = 1:numTypes
            count = count + 1;
            meanOnsetPert   = curSubj.audioMf0MeanPert_all{wC,ii}.ON.mean;
            meanOffsetPert  = curSubj.audioMf0MeanPert_all{wC,ii}.OF.mean;
            secPert(:,:,1) = meanOnsetPert;
            secPert(:,:,2) = meanOffsetPert;
            audioMf0SecPertAllSubj_all{wC,ii} = cat(2, audioMf0SecPertAllSubj_all{wC,ii}, secPert);

            meanOnsetPert   = curSubj.audioHf0MeanPert_all{wC,ii}.ON.mean;
            meanOffsetPert  = curSubj.audioHf0MeanPert_all{wC,ii}.OF.mean;
            secPert(:,:,1) = meanOnsetPert;
            secPert(:,:,2) = meanOffsetPert;
            audioHf0SecPertAllSubj_all{wC,ii} = cat(2, audioHf0SecPertAllSubj_all{wC,ii}, secPert);
        end
    end
end

allSubjRes.audioMf0SecCont = audioMf0SecContAllSubj;
allSubjRes.audioMf0SecPert = audioMf0SecPertAllSubj;
allSubjRes.audioHf0SecCont = audioHf0SecContAllSubj;
allSubjRes.audioHf0SecPert = audioHf0SecPertAllSubj;

allSubjRes.audioMf0SecPert_all = audioMf0SecPertAllSubj_all;
allSubjRes.audioMf0SecPert_all = audioHf0SecPertAllSubj_all;
end

function secDataM = meanSecData(secData)

OnsetSecs  = secData(:,:,1);
OffsetSecs = secData(:,:,2);
[~, numTrial] = size(OnsetSecs);

secDataM.ON.mean = nanmean(OnsetSecs, 2);
secDataM.OF.mean = nanmean(OffsetSecs, 2);

secDataM.ON.STD = nanstd(OnsetSecs, 0, 2);
secDataM.OF.STD = nanstd(OffsetSecs, 0, 2);

secDataM.ON.SEM = secDataM.ON.STD/sqrt(numTrial); % Standard Error
secDataM.OF.SEM = secDataM.OF.STD/sqrt(numTrial); % Standard Error

secDataM.ON.NCI = 1.96*secDataM.ON.SEM;  % 95% Confidence Interval
secDataM.OF.NCI = 1.96*secDataM.OF.SEM;  % 95% Confidence Interval
end

function audioDynamics_Somato = InflationResponse(secTime, secDataM)
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

onsetM       = secDataM.ON.mean;
[numSamp, ~] = size(onsetM); % Size of the data we are dealing with

ir = initInflationResponseStruct(); % Initialize the structure that handles the variable calculations
ir.time     = secTime;              % Time Interval for the sectioned trials (-0.5->1.0s)
ir.onset    = onsetM;   % f0 Trace sectioned around pert Onset.

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
ir.tAtRespRange = [0.8 1.0];              % Time Values in the period 800ms to 1000ms after perturbation onset
ir.iAtRespRange = find(ir.time >= ir.tAtRespRange(1) & ir.time<= ir.tAtRespRange(2));                % Last index in section
ir.vAtRespRange = ir.onset(ir.iAtRespRange);    % Time Value when participant 'fully responded' (1.0s)
ir.tAtResp      = mean(ir.tAtRespRange);
ir.vAtResp      = mean(ir.vAtRespRange);
ir.respMag = ir.vAtResp - ir.vAtMin; % Distance traveled from min f0 value to response f0 value

% RespPer
ir.respPer = 100*(ir.respMag/ir.stimMag); % Percent change from stimMag to respMag 

% Add to the audioDynamics struct
respVarM = [ir.tAtMin ir.stimMag ir.respMag ir.respPer];
audioDynamics_Somato.respVarM = respVarM;
% drawInflationResultMetrics(ir, 1, 0, 1)
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

function audioDynamics_Audio = PitchShiftReflexResponse(secTime, secDataM)
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

onsetM       = secDataM.ON.mean;
[numSamp, ~] = size(onsetM); % Size of the data we are dealing with

ir = initInflationResponseStruct(); % Initialize the structure that handles the variable calculations
ir.time     = secTime;              % Time Interval for the sectioned trials (-0.5->1.0s)
ir.onset    = onsetM;   % f0 Trace sectioned around pert Onset.

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
ir.tAtRespRange = [0.8 1.0];              % Time Values in the period 800ms to 1000ms after perturbation onset
ir.iAtRespRange = find(ir.time >= ir.tAtRespRange(1) & ir.time<= ir.tAtRespRange(2));                % Last index in section
ir.vAtRespRange = ir.onset(ir.iAtRespRange);    % Time Value when participant 'fully responded' (1.0s)
ir.tAtResp      = mean(ir.tAtRespRange);
ir.vAtResp      = mean(ir.vAtRespRange);
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

audioDynamics_Audio.tAtMin  = [];
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
    
    audioDynamics_Audio.tAtMin  = cat(1, audioDynamics_Audio.tAtMin, ir.tAtMin);
    audioDynamics_Audio.stimMag = cat(1, audioDynamics_Audio.stimMag, ir.stimMag);
    audioDynamics_Audio.respMag = cat(1, audioDynamics_Audio.respMag, ir.respMag);
    audioDynamics_Audio.respPer = cat(1, audioDynamics_Audio.respPer, ir.respPer);
end
audioDynamics_Audio.tAtMinM   = round(mean(audioDynamics_Audio.tAtMin),3);
audioDynamics_Audio.tAtMinSD  = round(mean(audioDynamics_Audio.tAtMin),3);
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

    [~, Imax] = max(audioMMean.ON.mean); %Max Pert Onset
    upBoundOn = round(audioMMean.ON.mean(Imax) + audioMMean.ON.NCI(Imax) + 10);
    [~, Imin] = min(audioMMean.ON.mean); %Min Pert Onset
    lwBoundOn = round(audioMMean.ON.mean(Imin) - audioMMean.ON.NCI(Imin) - 10);

    [~, Imax] = max(audioMMean.OF.mean); %Max Pert Offset
    upBoundOf = round(audioMMean.OF.mean(Imax) + audioMMean.OF.NCI(Imax) + 10);
    [~, Imin] = min(audioMMean.OF.mean); %Min Pert Offset
    lwBoundOf = round(audioMMean.OF.mean(Imin) - audioMMean.OF.NCI(Imin) - 10);

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
    
    [~, Imax] = max(audioHMean.ON.mean); %Max Pert Onset
    upBoundOn = round(audioHMean.ON.mean(Imax) + audioHMean.ON.NCI(Imax) + 10);
    [~, Imin] = min(audioHMean.ON.mean); %Min Pert Onset
    lwBoundOn = round(audioHMean.ON.mean(Imin) - audioHMean.ON.NCI(Imin) - 10);

    [~, Imax] = max(audioHMean.OF.mean); %Max Pert Offset
    upBoundOf = round(audioHMean.OF.mean(Imax) + audioHMean.OF.NCI(Imax) + 10);
    [~, Imin] = min(audioHMean.OF.mean); %Min Pert Offset
    lwBoundOf = round(audioHMean.OF.mean(Imin) - audioHMean.OF.NCI(Imin) - 10);

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
maxPres = max(ss.sensorPMean.ON.mean) + 0.5;
minPres = min(ss.sensorPMean.ON.mean) - 0.1;

lims.presMean = [-0.5 1.0 minPres maxPres];
end

function [adjustPres, InflDeflT] = scalePressureVals(ss)

time        = ss.secTimeP;
presOnsetM  = ss.sensorPMean.ON.mean;
presOnsetE  = ss.sensorPMean.ON.NCI;
presOffsetM = ss.sensorPMean.OF.mean;
presOffsetE = ss.sensorPMean.OF.NCI;

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

function polRes = packStatTableSingle(pA, polRes, sortStruc, curSubj)

for wC = 1:pA.numCond
    
    curIdx = (curSubj*pA.numCond) - pA.numCond + wC;
    polRes.statTableSingle.SubjID(curIdx)  = sortStruc.studyID;
    polRes.statTableSingle.Age(curIdx)     = sortStruc.age;
    polRes.statTableSingle.Gender(curIdx)  = sortStruc.gender;
    polRes.statTableSingle.f0(curIdx)      = sortStruc.f0b(wC);
    polRes.statTableSingle.AudFB(curIdx)   = sortStruc.AudFB{wC}(1);
    polRes.statTableSingle.StimMag(curIdx)   = sortStruc.respVarSingle{wC, 1}.stimMagM;
    polRes.statTableSingle.StimMagSD(curIdx) = sortStruc.respVarSingle{wC, 1}.stimMagSD;
    polRes.statTableSingle.RespMag(curIdx)   = sortStruc.respVarSingle{wC, 1}.respMagM;
    polRes.statTableSingle.RespMagSD(curIdx) = sortStruc.respVarSingle{wC, 1}.respMagSD;
    polRes.statTableSingle.RespPer(curIdx)   = sortStruc.respVarSingle{wC, 1}.respPerM;
    polRes.statTableSingle.RespPerSD(curIdx) = sortStruc.respVarSingle{wC, 1}.respPerSD;
end
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

function organizeAndSaveExcludedTrialTable(dirs, pA, allSubjRes, tossTrialTable, svFile)

fprintf('\nAcross all subjects, %d trials were excluded Automatically:\n', allSubjRes.tossedAll);
fprintf('%d trials due to late starts\n', allSubjRes.tossedLate);
fprintf('%d trials due to voice breaks\n', allSubjRes.tossedBreak);
fprintf('%d trials due to pitch miscalc\n', allSubjRes.tossedMisCalc);
fprintf('An additional %d trials were excluded Manually, of the %d trials Manually selected\n', allSubjRes.tossedAutoMiss, allSubjRes.tossedManual);
fprintf('Automatic trial exclusion methods accounted for %s%% of trials selected Manually\n', allSubjRes.autoSuccessPerc)
fprintf('\n')

if svFile == 1
    uf = uifigure('Position', [0 40 800 1200]);
    uitable(uf, 'Data', tossTrialTable{:,:},...
            'ColumnName', tossTrialTable.Properties.VariableNames,...
            'OuterPosition', [0 0 800 1200]);
        
    dirs.excludedTrialTable = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'ExcludedTrial.csv']);
    writetable(tossTrialTable, dirs.excludedTrialTable, 'WriteVariableNames',true)
end
end