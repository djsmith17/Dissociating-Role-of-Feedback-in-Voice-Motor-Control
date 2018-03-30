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
pA.pAnalysis     = 'LarynxPos'; % Change this name to load different pooled data sets

dirs               = dfDirs(pA.project);
dirs.SavResultsDir = fullfile(dirs.Results, 'Pooled Analyses', pA.pAnalysis);
dirs.PooledConfigF = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'PooledConfig.mat']);

% Can we find the Pooled Config File?
if exist(dirs.PooledConfigF, 'file') == 0
    fprintf('\nERROR: Pooled Config File %s does not exist! Please create it with a GenConfig Function\n', dirs.PooledConfigF)
    return
else
    % Load the configuration file. Should return a data structure cF
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

% allDataStr is 3D struc with dim (Parti nRun Cond);
allDataStr = [];
for ii = 1:pA.numPart
    participant = pA.participants{ii};

    subjRes  = [];
    condARes = []; condBRes = [];
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
            % Returns a results struture of res
        end
        
        % Which variable are we sorting against?
        condTest = eval(pA.condVar);
        % Which condition in our list, is the one in this run?
%         [~, condPos] = ismember(condTest, pA.cond);
        
        if strcmp(condTest, pA.cond{1})
            condARes = cat(2, condARes, res);
        else
            condBRes = cat(2, condBRes, res);
        end    
    end
    subjRes = cat(3, subjRes, condARes); % Cat the conditions along the z axis
    subjRes = cat(3, subjRes, condBRes); % Cat the conditions along the z axis
    
    allDataStr = cat(1, allDataStr, subjRes);
end

[~, numRunCond, ~] = size(allDataStr);

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
        
        thisStruc = initOrgStruct();
        thisStruc.subject = ['Participant ' num2str(ii)]; % Pooled Analysis Name
        thisStruc.curSess = [thisStruc.subject ' ' pA.cond{jj}];
        
        for kk = 1:numRunCond
            curRun = allDataStr(ii, kk, jj);
            
            thisStruc.studyID = curRun.subject; % Study ID
            thisStruc.AudFB   = curRun.AudFB;
            
            thisStruc.runs   = cat(1, thisStruc.runs, curRun.run);
            thisStruc.runf0b = cat(1, thisStruc.runf0b, curRun.f0b);
            
            thisStruc.allContTrials = cat(1, thisStruc.allContTrials, curRun.numContTrialsFin);
            thisStruc.allPertTrials = cat(1, thisStruc.allPertTrials, curRun.numPertTrialsFin);
            
            thisStruc.secTime = curRun.secTime;
            thisStruc.audioMf0SecPert = cat(2, thisStruc.audioMf0SecPert, curRun.audioMf0SecPert);
            thisStruc.audioMf0SecCont = cat(2, thisStruc.audioMf0SecCont, curRun.audioMf0SecCont);
            thisStruc.respVar         = cat(1, thisStruc.respVar, curRun.respVar);
        end
        
        thisStruc.f0b             = mean(thisStruc.runf0b);
        
        thisStruc.numContTrialsFin = sum(thisStruc.allContTrials);
        thisStruc.numPertTrialsFin = sum(thisStruc.allPertTrials);
        
        thisStruc.audioMf0MeanPert = meanRunAudioData(thisStruc.audioMf0SecPert);
        thisStruc.audioMf0MeanCont = meanRunAudioData(thisStruc.audioMf0SecCont);
        thisStruc.respVarM         = mean(thisStruc.respVar, 1);
        
        lims = identifyLimits(thisStruc, 0);
        thisStruc.limitsAmean = lims.audioMean;
        
        combDataStr(ii,jj) = thisStruc;        
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
allSubjRes.audioMf0MeanCont  = meanRunAudioData(allSubjRes.audioMf0SecCont);
allSubjRes.audioMf0MeanPertM = meanRunAudioData(allSubjRes.audioMf0SecPertM);
allSubjRes.audioMf0MeanPertV = meanRunAudioData(allSubjRes.audioMf0SecPertV);

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

dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'ResultsDRF.mat']);
fprintf('Saving Pooled Analysis for %s\n', pA.pAnalysis)
save(dirs.SavResultsFile, 'allDataStr', 'combDataStr', 'statLib', 'allSubjRes', 'statLibAll', 'pltNm')

dirs.excelFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'Stat.xlsx']);
xlswrite(dirs.excelFile, statLib, 1)
end

function thisStruc = initOrgStruct()

thisStruc.subject = [];
thisStruc.curSess = [];
thisStruc.studyID = [];
thisStruc.AudFB   = [];
thisStruc.runs    = {};

thisStruc.runf0b  = [];
thisStruc.f0b     = [];

thisStruc.allContTrials    = [];
thisStruc.numContTrialsFin = [];
thisStruc.allPertTrials    = [];
thisStruc.numPertTrialsFin = [];

thisStruc.secTime  = [];
thisStruc.audioMf0SecPert = [];
thisStruc.audioMf0SecCont = [];
thisStruc.respVar         = [];

thisStruc.audioMf0MeanPert = [];
thisStruc.audioMf0MeanCont = [];
thisStruc.respVarM         = [];

end

function meanAudio = meanRunAudioData(secAudio)

OnsetSecs  = secAudio(:,:,1);
OffsetSecs = secAudio(:,:,2);
[~, numTrial] = size(OnsetSecs);

meanOnset  = mean(OnsetSecs, 2);
meanOffset = mean(OffsetSecs, 2);

stdOnset   = std(OnsetSecs, 0, 2);
stdOffset  = std(OffsetSecs, 0, 2);

SEMOnset   = stdOnset/sqrt(numTrial);  % Standard Error
SEMOffset  = stdOffset/sqrt(numTrial); % Standard Error

NCIOnset   = 1.96*SEMOnset;  % 95% Confidence Interval
NCIOffset  = 1.96*SEMOffset; % 95% Confidence Interval

meanAudio = [meanOnset NCIOnset meanOffset NCIOffset];
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

%Full Inidividual Trials: Pressure Sensor
lims.pressure   = [0 4 0 5];

%Aligned Pressure Data
lims.pressureAl = [0 3.5 0 5];

%Full Individual Trials: Force Sensors
lims.force      = [0 4 1 5];

%Full trial f0 analysis
%Full Individual Trials: f0 Audio 
lims.audio      = [0 4 -100 100];

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