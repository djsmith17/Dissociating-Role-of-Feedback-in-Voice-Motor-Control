function dfRunPooledAnalysis()
%Analyses the microphone data from the somatosensory perturbation
%experiment. Measures the change in f0 over each trial, and each run for a
%given participant. At the end it approximates a general response to
%inflation to be used in the auditory perturbation experiment

%Require the Signal Processing Toolbox

clear all; close all; clc
pA.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control'; 
pA.pAnalysis     = 'SfN2017';
pA.participants  = {'Pilot24'; 'Pilot25'; 'Pilot26'; 'Pilot22'}; %List of multiple participants.
pA.numPart       = length(pA.participants);
pA.runs          = {'SF1'; 'SF2'; 'SF3'; 'SF4'}; %All runs to consider 
pA.numRuns       = length(pA.runs);
pA.cond          = {' Masking Noise'; ' Normal Voicing'};

dirs               = dfDirs(pA.project);
dirs.SavResultsDir = fullfile(dirs.Results, 'Pooled Analyses', pA.pAnalysis);

if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

% allDataStr is 3D struc with dim (Parti nRun Cond);
for ii = 1:pA.numPart
    participant = pA.participants{ii};
    colM = 1; colV = 1;
    fprintf('Sorting Runs for %s\n', participant)
    for jj = 1:pA.numRuns
        if ii == 1 & jj == 3
            disp('lol')
        else
            run         = pA.runs{jj};
            dirs.SavFileDir  = fullfile(dirs.Results, participant, run); %Where to save results        
            dirs.SavFile     = fullfile(dirs.SavFileDir, [participant run 'ResultsDRF.mat']);

            if exist(dirs.SavFile, 'file') == 0
                disp('ERROR: NO DANG FILE')
                return
            end        
            load(dirs.SavFile)
            pA.AudFB = niRes.AudFB;

            if strcmp(pA.AudFB, 'Masking Noise')
               allDataStr(ii, colM, 1) = niRes;
               colM = colM + 1;
            else
               allDataStr(ii, colV, 2) = niRes;
               colV = colV + 1;
            end    
        end
    end
end

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
    fprintf('Combining task conditions for %s\n',participant)
    for jj = 1:2 %Masking Noise, then Voice Conditions
        runSt1 = allDataStr(ii, 1, jj);
        runSt2 = allDataStr(ii, 2, jj);
        
        thisStruc.parti           = runSt1.subject;
        thisStruc.subject         = ['Participant ' num2str(ii)]; %doubleblind sorta, I guess. Shoot me
        thisStruc.runs            = {runSt1.run; runSt2.run};
        thisStruc.curSess         = [thisStruc.subject pA.cond{jj}];
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