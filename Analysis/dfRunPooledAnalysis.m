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
pA.pAnalysis     = 'SfN2017'; % Change this name to load different pooled data sets Ex: SfN2017, LarynxPos

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
[~, pA.numRuns]  = size(pA.runs);
pA.cond          = cF.cond;         % Conditions to test against
pA.numCond       = length(pA.cond); 
pA.condVar       = cF.condVar;      % Variable to test the condition

pA.pltNameMVi    = cF.pltNameMVi;
pA.pltNameMVm    = cF.pltNameMVm;

% Load all saved results and order into a large data structure
allDataStr = []; % (numPart x numRuns)
for ii = 1:pA.numPart
    participant = pA.participants{ii};
    fprintf('Loading Runs for %s\n', participant)
    
    subjRes  = [];
    for jj = 1:pA.numRuns
        run              = pA.runs{ii, jj};
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
        
        sortStruc  = combineCondTrials(pA, curRes, sortStruc);       
        allSubjRes = combineCondTrials(pA, curRes, allSubjRes);
    end
        
    sortStruc = meanCondTrials(pA, sortStruc);
    sortStruc.pltName = pA.pltNameMVi{ii};
    
    pooledRunStr(ii)   = sortStruc;        
end

allSubjRes = meanCondTrials(pA, allSubjRes);
allSubjRes.pltName  = pA.pltNameMVm;

% Save the Pooled Results
dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'ResultsDRF.mat']);
fprintf('Saving Pooled Analysis for %s\n', pA.pAnalysis)
save(dirs.SavResultsFile, 'pooledRunStr', 'allSubjRes')

dirs.excelFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'Stat.xlsx']);
% xlswrite(dirs.excelFile, statLib, 1)
end

function sortStr = initSortedStruct(numCond)

sortStr.subject = [];
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
sortStr.respVar         = cell(numCond, 1);

sortStr.audioMf0MeanPert = cell(numCond, 1);
sortStr.audioMf0MeanCont = [];
sortStr.respVarM         = zeros(numCond, 4);

sortStr.tossedAll        = [];
sortStr.tossedLate       = [];
sortStr.tossedBreak      = [];
sortStr.tossedMisCalc    = [];
end

function polRes = combineCondTrials(pA, curRes, polRes)

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
end

function polRes = meanCondTrials(pA, polRes)

polRes.numContTrialsFin = sum(polRes.allContTrials);
polRes.audioMf0MeanCont = meanSecData(polRes.audioMf0SecCont);
for kk = 1:pA.numCond
    polRes.f0b(kk)              = mean(polRes.runf0b{kk});
    
    polRes.numPertTrialsFin(kk) = sum(polRes.allPertTrials{kk});
    polRes.audioMf0MeanPert{kk} = meanSecData(polRes.audioMf0SecPert{kk});
    polRes.respVarM(kk, :)      = mean(polRes.respVar{kk}, 1);
end

lims = identifyLimits(polRes);
polRes.limitsAmean = lims.audioMean;

statLib         = packStatLib(polRes);
polRes.statLib  = statLib;
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

function lims = identifyLimits(ss)

mf0MeanPert = ss.audioMf0MeanPert;
numCond = length(mf0MeanPert);

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

maxUpBound = max(setupBoundSec);
minLwBound = min(setlwBoundSec);

lims.audioMean = [-0.5 1.0 minLwBound maxUpBound];
end