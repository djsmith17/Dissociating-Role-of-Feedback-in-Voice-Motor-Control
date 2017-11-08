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
        run         = pA.runs{jj};
        dirs.SavFileDir  = fullfile(dirs.Results, participant, run); %Where to save results        
        dirs.SavFile     = fullfile(dirs.SavFileDir, [participant run 'ResultsDRF.mat']);
        
        if exist(dirs.SavFile, 'file') == 0
            disp('ERROR: NO DANG FILE')
            return
        end        
        load(dirs.SavFile)
        AudFB = auAn.AudFB;
        
        if strcmp(AudFB, 'Masking Noise')
           allDataStr(ii, colM, 1) = niRes;
           colM = colM + 1;
        else
           allDataStr(ii, colV, 2) = niRes;
           colV = colV + 1;
        end       
    end
end

cond = {' Masking Noise'; ' Normal Voicing'};
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
        thisStruc.curSess         = [thisStruc.subject cond{jj}];
        thisStruc.numContTrials   = sum([runSt1.numContTrials runSt2.numContTrials]);
        thisStruc.numPertTrials   = sum([runSt1.numPertTrials runSt2.numPertTrials]);
        thisStruc.secTime         = runSt1.secTime;
        thisStruc.runf0b          = [runSt1.f0b runSt2.f0b];
        thisStruc.audioMf0SecPert = [runSt1.audioMf0SecPert runSt2.audioMf0SecPert];
        thisStruc.audioMf0SecCont = [runSt1.audioMf0SecCont runSt2.audioMf0SecCont];        
        thisStruc.respVar         = [runSt1.respVar; runSt2.respVar];
        
        thisStruc.f0b              = mean(thisStruc.runf0b);
        thisStruc.audioMf0MeanPert = meanRunAudioData(thisStruc.audioMf0SecPert);
        thisStruc.audioMf0MeanCont = meanRunAudioData(thisStruc.audioMf0SecCont);
        thisStruc.respVarm         = mean(thisStruc.respVar, 1);
        
        lims = identifyLimits(thisStruc);
        thisStruc.limitsAmean = lims.audioMean;
        
        combDataStr(ii,jj) = thisStruc;        
    end
    mask = combDataStr(ii,1);
    voic = combDataStr(ii,2); 
    
    if mask.numPertTrials ~= voic.numPertTrials
        if mask.numPertTrials > voic.numPertTrials
            disp('There are less Voiced Trials')
            t2u = 1:voic.numPertTrials;
        else
            disp('There are less Masked Trials')
            t2u = 1:mask.numPertTrials; 
        end
    else
        disp('Using all pert trials')
        t2u = 1:mask.numPertTrials;
    end
    
    [Hstim, pStim] = ttest(mask.respVar(t2u,2), voic.respVar(t2u,2));
    [Hresp, pResp] = ttest(mask.respVar(t2u,3), voic.respVar(t2u,3));
    [Hperc, pPerc] = ttest(mask.respVar(t2u,4), voic.respVar(t2u,4));
    
    statLib(ii,1) = mask.respVarm(2); %Masking StimMag
    statLib(ii,2) = voic.respVarm(2); %Voicing StimMag
    statLib(ii,3) = mask.respVarm(3); %Masking RespMag
    statLib(ii,4) = voic.respVarm(3); %Voicing RespMag
    statLib(ii,5) = mask.respVarm(4); %Masking %
    statLib(ii,6) = voic.respVarm(4); %Voicing %
    statLib(ii,7) = pStim; %p-value stimulus
    statLib(ii,8) = pResp; %p-value response
    statLib(ii,9) = pPerc; %p-value percent increase 
end

% statTable = table(statLib(:,1), 

dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'ResultsDRF.mat']);
fprintf('Saving Pooled Analysis for %s\n', pA.pAnalysis)
save(dirs.SavResultsFile, 'allDataStr', 'combDataStr', 'statLib')

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

function lims = identifyLimits(niAn)

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
[~, Imax] = max(niAn.audioMf0MeanPert(:,1)); %Max Pert Onset
upBoundOn = round(niAn.audioMf0MeanPert(Imax,1) + niAn.audioMf0MeanPert(Imax,2) + 10);
[~, Imin] = min(niAn.audioMf0MeanPert(:,1)); %Min Pert Onset
lwBoundOn = round(niAn.audioMf0MeanPert(Imin,1) - niAn.audioMf0MeanPert(Imin,2) - 10);

[~, Imax] = max(niAn.audioMf0MeanPert(:,3)); %Max Pert Offset
upBoundOf = round(niAn.audioMf0MeanPert(Imax,3) + niAn.audioMf0MeanPert(Imax,4) + 10);
[~, Imin] = min(niAn.audioMf0MeanPert(:,3)); %Min Pert Offset
lwBoundOf = round(niAn.audioMf0MeanPert(Imin,3) - niAn.audioMf0MeanPert(Imin,4) - 10);

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