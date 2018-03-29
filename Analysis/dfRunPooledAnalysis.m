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
pA.pAnalysis     = 'SfN2017'; % Change this name to load different pooled data sets

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
            % Returns a results struture of niRes
        end
        
        % Which variable are we sorting against?
        condTest = eval(pA.condVar);
        % Which condition in our list, is the one in this run?
%         [~, condPos] = ismember(condTest, pA.cond);
        
        if strcmp(condTest, pA.cond{1})
            condARes = cat(2, condARes, niRes);
        else
            condBRes = cat(2, condBRes, niRes);
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
allSubjRes.audioMf0SecPertM = [];
allSubjRes.audioMf0SecPertV = [];
allSubjRes.audioMf0SecContM = [];
allSubjRes.audioMf0SecContV = [];
unSubM.respVar           = [];
unSubV.respVar           = [];

statLib = [];
for ii = 1:pA.numPart
    participant = pA.participants{ii};
    fprintf('Combining task conditions for %s\n', participant)
    for jj = 1:2 % Masking Noise, then Voice Conditions
        runSt1 = allDataStr(ii, 1, jj);
        runSt2 = allDataStr(ii, 2, jj);
        
        thisStruc.studyID         = runSt1.subject;               % Study ID
        thisStruc.subject         = ['Participant ' num2str(ii)]; % Pooled Analysis Name
        thisStruc.runs            = {runSt1.run; runSt2.run};
        thisStruc.curSess         = [thisStruc.subject ' ' pA.cond{jj}];
        thisStruc.AudFB           = runSt1.AudFB;
        thisStruc.numContTrials   = sum([runSt1.numContTrials runSt2.numContTrialsPP]);
        thisStruc.numPertTrials   = sum([runSt1.numPertTrials runSt2.numPertTrialsPP]);
        
        thisStruc.secTime         = runSt1.secTime;
        thisStruc.runf0b          = [runSt1.f0b runSt2.f0b];
        thisStruc.audioMf0SecPert = [runSt1.audioMf0SecPert runSt2.audioMf0SecPert];
        thisStruc.audioMf0SecCont = [runSt1.audioMf0SecCont runSt2.audioMf0SecCont];        
        thisStruc.respVar         = [runSt1.respVar; runSt2.respVar];
        
        thisStruc.f0b              = mean(thisStruc.runf0b);
        thisStruc.audioMf0MeanPert = meanRunAudioData(thisStruc.audioMf0SecPert);
        thisStruc.audioMf0MeanCont = meanRunAudioData(thisStruc.audioMf0SecCont);
        thisStruc.respVarm         = mean(thisStruc.respVar, 1);
        
        lims = identifyLimits(thisStruc, 0);
        thisStruc.limitsAmean = lims.audioMean;
        
        combDataStr(ii,jj) = thisStruc;        
    end
    mask = combDataStr(ii,1);
    voic = combDataStr(ii,2); 
    
    statLib(ii,:) = packStatLib(mask, voic);
    
    allSubjRes.numControlTrials = allSubjRes.numControlTrials + mask.numContTrials + voic.numContTrials;
    allSubjRes.numMaskedTrials = allSubjRes.numMaskedTrials + mask.numPertTrials;
    allSubjRes.numVoicedTrials = allSubjRes.numVoicedTrials + voic.numPertTrials;
    
    allSubjRes.audioMf0SecPertM = cat(2, allSubjRes.audioMf0SecPertM, mask.audioMf0SecPert);
    allSubjRes.audioMf0SecPertV = cat(2, allSubjRes.audioMf0SecPertV, voic.audioMf0SecPert);
    
    allSubjRes.audioMf0SecContM = cat(2, allSubjRes.audioMf0SecContM, mask.audioMf0SecCont);
    allSubjRes.audioMf0SecContV = cat(2, allSubjRes.audioMf0SecContV, voic.audioMf0SecCont);
    
    unSubM.respVar = cat(1, unSubM.respVar, mask.respVar);
    unSubV.respVar = cat(1, unSubV.respVar, voic.respVar);
end
allSubjRes.secTime           = mask.secTime;
allSubjRes.audioMf0MeanPertM = meanRunAudioData(allSubjRes.audioMf0SecPertM);
allSubjRes.audioMf0MeanPertV = meanRunAudioData(allSubjRes.audioMf0SecPertV);
allSubjRes.audioMf0MeanContM = meanRunAudioData(allSubjRes.audioMf0SecContM);
allSubjRes.audioMf0MeanContV = meanRunAudioData(allSubjRes.audioMf0SecContV);

unSubM.respVarm              = mean(unSubM.respVar, 1);
unSubV.respVarm              = mean(unSubV.respVar, 1);

allSubjRes.respVarM          = unSubM.respVar;
allSubjRes.respVarV          = unSubV.respVar;
allSubjRes.respVarmM         = unSubM.respVarm;
allSubjRes.respVarmV         = unSubV.respVarm;

limsM = identifyLimits(allSubjRes, 1);
allSubjRes.limitsAmeanM = limsM.audioMean;
limsV = identifyLimits(allSubjRes, 2);
allSubjRes.limitsAmeanV = limsV.audioMean;

statLibAll = packStatLib(unSubM, unSubV);

dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'ResultsDRF.mat']);
fprintf('Saving Pooled Analysis for %s\n', pA.pAnalysis)
save(dirs.SavResultsFile, 'allDataStr', 'combDataStr', 'statLib', 'allSubjRes', 'statLibAll')

dirs.excelFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'Stat.xlsx']);
xlswrite(dirs.excelFile, statLib, 1)
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

statLib(1) = mask.respVarm(2); %Masking StimMag
statLib(2) = voic.respVarm(2); %Voicing StimMag
statLib(3) = mask.respVarm(3); %Masking RespMag
statLib(4) = voic.respVarm(3); %Voicing RespMag
statLib(5) = mask.respVarm(4); %Masking %
statLib(6) = voic.respVarm(4); %Voicing %
statLib(7) = pStim; %p-value stimulus
statLib(8) = pResp; %p-value response
statLib(9) = pPerc; %p-value percent increase 
end

function lims = identifyLimits(niAn, fl)

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
    audioMean = niAn.audioMf0MeanPertM;
elseif fl == 2
    audioMean = niAn.audioMf0MeanPertV;
else
    audioMean = niAn.audioMf0MeanPert;
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