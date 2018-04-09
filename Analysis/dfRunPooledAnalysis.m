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

pltNm.pltNameMVi  = cF.pltNameMVi;
pltNm.pltNameMVm  = cF.pltNameMVm;

% Load all saved results and order into a large data structure
allDataStr = []; % (numPart x numRuns)
for ii = 1:pA.numPart
    participant = pA.participants{ii};

    subjRes  = [];
    fprintf('Sorting Runs for %s\n', participant)
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

allSubjRes.numControlTrials = 0;
allSubjRes.numMaskedTrials  = 0;
allSubjRes.numVoicedTrials  = 0;
allSubjRes.secTime          = [];
allSubjRes.audioMf0SecCont  = [];
allSubjRes.audioMf0SecPertM = [];
allSubjRes.audioMf0SecPertV = [];
unSubM.respVar           = [];
unSubV.respVar           = [];

statLib = [];
for ii = 1:pA.numPart
    participant = pA.participants{ii};
    fprintf('Combining task conditions for %s\n', participant)
    
    
    for jj = 1:2 % Masking Noise, then Voice Conditions
        
        unpkStruc = initUnpackStruct();
        unpkStruc.subject = ['Participant ' num2str(ii)]; % Pooled Analysis Name
        unpkStruc.curSess = [unpkStruc.subject ' ' pA.cond{jj}];
        
        for kk = 1:numRunCond
            curRun = allDataStr(ii, kk, jj);
            
            unpkStruc.studyID = curRun.subject; % Study ID
            unpkStruc.AudFB   = curRun.AudFB;
            
            unpkStruc.runs   = cat(1, unpkStruc.runs, curRun.run);
            unpkStruc.runf0b = cat(1, unpkStruc.runf0b, curRun.f0b);
            
            unpkStruc.allContTrials = cat(1, unpkStruc.allContTrials, curRun.numContTrialsFin);
            unpkStruc.allPertTrials = cat(1, unpkStruc.allPertTrials, curRun.numPertTrialsFin);
            
            unpkStruc.secTime = curRun.secTime;
            unpkStruc.audioMf0SecPert = cat(2, unpkStruc.audioMf0SecPert, curRun.audioMf0SecPert);
            unpkStruc.audioMf0SecCont = cat(2, unpkStruc.audioMf0SecCont, curRun.audioMf0SecCont);
            unpkStruc.respVar         = cat(1, unpkStruc.respVar, curRun.respVar);
        end
        
        unpkStruc.f0b             = mean(unpkStruc.runf0b);
        
        unpkStruc.numContTrialsFin = sum(unpkStruc.allContTrials);
        unpkStruc.numPertTrialsFin = sum(unpkStruc.allPertTrials);
        
        unpkStruc.audioMf0MeanPert = meanSecData(unpkStruc.audioMf0SecPert);
        unpkStruc.audioMf0MeanCont = meanSecData(unpkStruc.audioMf0SecCont);
        unpkStruc.respVarM         = mean(unpkStruc.respVar, 1);
        
        lims = identifyLimits(unpkStruc, 0);
        unpkStruc.limitsAmean = lims.audioMean;
        
        combDataStr(ii,jj) = unpkStruc;        
    end
    mask = combDataStr(ii,1);
    voic = combDataStr(ii,2); 
    
    statLib(ii,:) = packStatLib(mask, voic);
    
    allSubjRes.numControlTrials = allSubjRes.numControlTrials + mask.numContTrialsFin + voic.numContTrialsFin;
    allSubjRes.numMaskedTrials = allSubjRes.numMaskedTrials + mask.numPertTrialsFin;
    allSubjRes.numVoicedTrials = allSubjRes.numVoicedTrials + voic.numPertTrialsFin;
    
    allSubjRes.audioMf0SecPertM = cat(2, allSubjRes.audioMf0SecPertM, mask.audioMf0SecPert);
    allSubjRes.audioMf0SecPertV = cat(2, allSubjRes.audioMf0SecPertV, voic.audioMf0SecPert);
    
    % This will take all the control trials from all conditions and
    % concatenate them in one big matrix
    allSubjRes.audioMf0SecCont = cat(2, allSubjRes.audioMf0SecCont, mask.audioMf0SecCont);
    allSubjRes.audioMf0SecCont = cat(2, allSubjRes.audioMf0SecCont, voic.audioMf0SecCont);
    
    unSubM.respVar = cat(1, unSubM.respVar, mask.respVar);
    unSubV.respVar = cat(1, unSubV.respVar, voic.respVar);
end

allSubjRes.secTime           = mask.secTime;
allSubjRes.audioMf0MeanCont  = meanSecData(allSubjRes.audioMf0SecCont);
allSubjRes.audioMf0MeanPertM = meanSecData(allSubjRes.audioMf0SecPertM);
allSubjRes.audioMf0MeanPertV = meanSecData(allSubjRes.audioMf0SecPertV);

unSubM.respVarM              = mean(unSubM.respVar, 1);
unSubV.respVarM              = mean(unSubV.respVar, 1);

allSubjRes.respVarM          = unSubM.respVar;
allSubjRes.respVarV          = unSubV.respVar;
allSubjRes.respVarmM         = unSubM.respVarM;
allSubjRes.respVarmV         = unSubV.respVarM;

limsM = identifyLimits(allSubjRes, 1);
allSubjRes.limitsAmeanM = limsM.audioMean;
limsV = identifyLimits(allSubjRes, 2);
allSubjRes.limitsAmeanV = limsV.audioMean;

statLibAll = packStatLib(unSubM, unSubV);

if strcmp(pA.pAnalysis, 'LarynxPos') == 1
    [CRi, CRm] = collarResultConcat(allDataStr);
else
    CRi = []; CRm = [];
end


% Save the Pooled Results
dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'ResultsDRF.mat']);
fprintf('Saving Pooled Analysis for %s\n', pA.pAnalysis)
save(dirs.SavResultsFile, 'allDataStr', 'combDataStr', 'statLib', 'allSubjRes', 'statLibAll', 'pltNm', 'CRi', 'CRm')

dirs.excelFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'Stat.xlsx']);
% xlswrite(dirs.excelFile, statLib, 1)
end

function unpkStr = initUnpackStruct()

unpkStr.subject = [];
unpkStr.curSess = [];
unpkStr.studyID = [];
unpkStr.AudFB   = [];
unpkStr.runs    = {};

unpkStr.runf0b  = [];
unpkStr.f0b     = [];

unpkStr.allContTrials    = [];
unpkStr.numContTrialsFin = [];
unpkStr.allPertTrials    = [];
unpkStr.numPertTrialsFin = [];

unpkStr.secTime         = [];
unpkStr.audioMf0SecPert = [];
unpkStr.audioMf0SecCont = [];
unpkStr.respVar         = [];

unpkStr.audioMf0MeanPert = [];
unpkStr.audioMf0MeanCont = [];
unpkStr.respVarM         = [];

unpkStr.tossedAll        = [];
unpkStr.tossedLate       = [];
unpkStr.tossedBreak      = [];
unpkStr.tossedMisCalc    = [];
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

function [CRi, CRm] = collarResultConcat(allDataStr)

[numSubj, numCollarPos, numCond] = size(allDataStr);
nameCollarPos = {'LP', 'uLP', 'CC'};

CRm.curSess          = 'MeanSubjResp_CollarLoc';
CRm.numContTrials    = 0;
CRm.numPertTrialsLP  = 0;
CRm.numPertTrialsuLP = 0;
CRm.numPertTrialsCC  = 0;

CRm.secTime          = [];

CRm.audioMf0SecCont    = [];
CRm.audioMf0SecPertLP  = [];
CRm.audioMf0SecPertuLP = [];
CRm.audioMf0SecPertCC  = [];
unSubLP.respVar        = [];
unSubuLP.respVar       = [];
unSubCC.respVar        = [];

for ii = 1:numSubj
    for kk = 1:numCollarPos
        subjColl = initUnpackStruct();
        
        subjColl.subject = allDataStr(ii, 1, 1).subject;
        subjColl.curSess = [subjColl.subject 'Resp_CollarLoc'];
        subjColl.CollarPos = nameCollarPos{kk};
    
        for jj = 1:numCond %Masked then Voiced
            curRun = allDataStr(ii, kk, jj);

            subjColl.runs            = cat(1, subjColl.runs, curRun.run);
            subjColl.runf0b          = cat(1, subjColl.runf0b, curRun.f0b);
            
            subjColl.allContTrials = cat(1, subjColl.allContTrials, curRun.numContTrialsFin);
            subjColl.allPertTrials = cat(1, subjColl.allPertTrials, curRun.numPertTrialsFin);
            
            subjColl.secTime = curRun.secTime;
            subjColl.audioMf0SecPert = cat(2, subjColl.audioMf0SecPert, curRun.audioMf0SecPert);
            subjColl.audioMf0SecCont = cat(2, subjColl.audioMf0SecCont, curRun.audioMf0SecCont);
            subjColl.respVar         = cat(1, subjColl.respVar, curRun.respVar);
            
            numTossed                = length(curRun.removedTrialTracker);
            subjColl.tossedAll       = cat(1, subjColl.tossedAll, numTossed);
            
            if numTossed > 0  
                idxLateC = strfind(curRun.removedTrialTracker(:,2), 'Participant started too late!!');
                [~, idxLate] = find(not(cellfun('isempty', idxLateC)));
                
                idxBreakC = strfind(curRun.removedTrialTracker(:,2), 'Participant had a voice break!!');
                [~, idxBreak] = find(not(cellfun('isempty', idxBreakC)));
                
                idxMisC = strfind(curRun.removedTrialTracker(:,2), 'Miscalculated pitch Trace');
                [~, idxMis] = find(not(cellfun('isempty', idxMisC)));   
                
                numTossedLate  = sum(idxLate);
                numTossedBreak = sum(idxBreak);
                numTossedMisCalc  = sum(idxMis);       
            else
                numTossedLate = 0;
                numTossedBreak = 0;
                numTossedMisCalc  = 0;
            end
            subjColl.tossedLate       = cat(1, subjColl.tossedLate, numTossedLate);
            subjColl.tossedBreak      = cat(1, subjColl.tossedBreak, numTossedBreak);
            subjColl.tossedMisCalc    = cat(1, subjColl.tossedMisCalc, numTossedMisCalc);      
        end
        
        subjColl.f0b             = mean(subjColl.runf0b);
        
        subjColl.numContTrialsFin = sum(subjColl.allContTrials);
        subjColl.numPertTrialsFin = sum(subjColl.allPertTrials);
        
        subjColl.audioMf0MeanPert = meanSecData(subjColl.audioMf0SecPert);
        subjColl.audioMf0MeanCont = meanSecData(subjColl.audioMf0SecCont);
        subjColl.respVarM         = mean(subjColl.respVar, 1);
        
        subjColl.perTossed        = round(100*(sum(subjColl.tossedAll)/20), 1);
        subjColl.perTossedLate    = round(100*(sum(subjColl.tossedLate)/20), 1);
        subjColl.perTossedBreak   = round(100*(sum(subjColl.tossedBreak)/20), 1);
        subjColl.perTossedMisCalc = round(100*(sum(subjColl.tossedMisCalc)/20), 1);
        
        lims = identifyLimits(subjColl, 0);
        subjColl.limitsAmean = lims.audioMean;     
        
        collectedSet(ii, kk) = subjColl;        
    end
    
    LP  = collectedSet(ii,1);
    uLP = collectedSet(ii,2); 
    CC  = collectedSet(ii,3);
    
    CRi(ii).curSess          = LP.curSess;
    CRi(ii).numContTrials    = LP.numContTrialsFin + uLP.numContTrialsFin + CC.numContTrialsFin;
    CRi(ii).numPertTrialsLP  = LP.numPertTrialsFin;
    CRi(ii).numPertTrialsuLP = uLP.numPertTrialsFin;
    CRi(ii).numPertTrialsCC  = CC.numPertTrialsFin;

    CRi(ii).secTime             = LP.secTime;
    CRi(ii).audioMf0MeanCont    = meanSecData([LP.audioMf0SecCont, uLP.audioMf0SecCont, CC.audioMf0SecCont]);
    CRi(ii).audioMf0MeanPertLP  = LP.audioMf0MeanPert;
    CRi(ii).audioMf0MeanPertuLP = uLP.audioMf0MeanPert;
    CRi(ii).audioMf0MeanPertCC  = CC.audioMf0MeanPert;
    
%     identifyLimits = 
%     CRi(ii).limitsAMean  = 

%     statLib(ii,:) = packStatLib(mask, voic);
%     
    CRm.numContTrials    = CRm.numContTrials + LP.numContTrialsFin + uLP.numContTrialsFin + CC.numContTrialsFin;
    CRm.numPertTrialsLP  = CRm.numPertTrialsLP + LP.numPertTrialsFin;
    CRm.numPertTrialsuLP = CRm.numPertTrialsuLP + uLP.numPertTrialsFin;
    CRm.numPertTrialsCC  = CRm.numPertTrialsCC + CC.numPertTrialsFin;
    
    CRm.audioMf0SecPertLP  = cat(2, CRm.audioMf0SecPertLP, LP.audioMf0SecPert);
    CRm.audioMf0SecPertuLP = cat(2, CRm.audioMf0SecPertuLP, uLP.audioMf0SecPert);
    CRm.audioMf0SecPertCC  = cat(2, CRm.audioMf0SecPertCC, CC.audioMf0SecPert);
    
    % This will take all the control trials from all conditions and
    % concatenate them in one big matrix
    CRm.audioMf0SecCont = cat(2, CRm.audioMf0SecCont, LP.audioMf0SecCont);
    CRm.audioMf0SecCont = cat(2, CRm.audioMf0SecCont, uLP.audioMf0SecCont);
    CRm.audioMf0SecCont = cat(2, CRm.audioMf0SecCont, CC.audioMf0SecCont);
    
    unSubLP.respVar  = cat(1, unSubLP.respVar, LP.respVar);
    unSubuLP.respVar = cat(1, unSubuLP.respVar, uLP.respVar);
    unSubCC.respVar  = cat(1, unSubCC.respVar, CC.respVar);
end

CRm.secTime             = LP.secTime;
CRm.audioMf0MeanCont    = meanSecData(CRm.audioMf0SecCont);
CRm.audioMf0MeanPertLP  = meanSecData(CRm.audioMf0SecPertLP);
CRm.audioMf0MeanPertuLP = meanSecData(CRm.audioMf0SecPertuLP);
CRm.audioMf0MeanPertCC  = meanSecData(CRm.audioMf0SecPertCC);

unSubLP.respVarM      = mean(unSubLP.respVar, 1);
unSubuLP.respVarM     = mean(unSubuLP.respVar, 1);
unSubCC.respVarM      = mean(unSubCC.respVar, 1);

CRm.respVarLP          = unSubLP.respVar;
CRm.respVaruLP         = unSubuLP.respVar;
CRm.respVarCC          = unSubCC.respVar;
CRm.respVarMLP         = unSubLP.respVarM;
CRm.respVarMuLP        = unSubuLP.respVarM;
CRm.respVarMCC         = unSubCC.respVarM;
% 
% limsLP = identifyLimits(CRm, 1);
% CRm.limitsAmeanM = limsLP.audioMean;
% limsuLP = identifyLimits(CRm, 2);
% CRm.limitsAmeanV = limsuLP.audioMean;
% limsCC = identifyLimits(CRm, 2);
% CRm.limitsAmeanV = limsCC.audioMean;
end

function statLib = packStatLib(mask, voic)

[~, pStim] = ttest2(mask.respVar(:,2), voic.respVar(:,2));
[~, pResp] = ttest2(mask.respVar(:,3), voic.respVar(:,3));
[~, pPerc] = ttest2(mask.respVar(:,4), voic.respVar(:,4));

statLib(1) = mask.respVarM(2); %Masking StimMag
statLib(2) = voic.respVarM(2); %Voicing StimMag
statLib(3) = mask.respVarM(3); %Masking RespMag
statLib(4) = voic.respVarM(3); %Voicing RespMag
statLib(5) = mask.respVarM(4); %Masking %
statLib(6) = voic.respVarM(4); %Voicing %
statLib(7) = pStim; %p-value stimulus
statLib(8) = pResp; %p-value response
statLib(9) = pPerc; %p-value percent increase 
end

function lims = identifyLimits(An, fl)

%Section Mean Pertrubed Trials: f0 Audio 
if fl == 1
    audioMean = An.audioMf0MeanPertM;
elseif fl == 2
    audioMean = An.audioMf0MeanPertV;
else
    audioMean = An.audioMf0MeanPert;
end

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

lims.audioMean = [-0.5 1.0 lwBoundSec upBoundSec];
end