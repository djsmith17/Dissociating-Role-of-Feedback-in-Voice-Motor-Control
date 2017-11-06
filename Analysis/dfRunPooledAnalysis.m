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
    colM = 1; colV = 1;
    for jj = 1:pA.numRuns
        participant = pA.participants{ii};
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
    for jj = 1:2 %Masking Noise, then Voice Conditions
        runSt1 = allDataStr(ii, 1, jj);
        runSt2 = allDataStr(ii, 2, jj);
        
        thisStruc.subject         = runSt1.subject;
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
        thisStruc.limitsAmean      = [-0.5 1.0 -80 20];
        thisStruc.respVarm         = mean(thisStruc.respVar, 1);
        
        combDataStr(ii,jj) = thisStruc;        
    end
    mask = combDataStr(ii,1);
    voic = combDataStr(ii,2); 
    
    [Hstim, pStim] = ttest(mask.respVar(:,2), voic.respVar(:,2));
    [Hresp, pResp] = ttest(mask.respVar(:,3), voic.respVar(:,3));
    [Hperc, pPerc] = ttest(mask.respVar(:,4), voic.respVar(:,4));
    
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

dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'ResultsDRF.mat']);
fprintf('Saving Pooled Analysis for %s\n', pA.pAnalysis)
save(dirs.SavResultsFile, 'allDataStr', 'combDataStr', 'statLib')
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