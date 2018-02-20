function An = dfAnalysisAudio(dirs, An, AudFlag, varargin)
%I found that I was mostly writing the same things twice in my scripts for
%audapter analysis and NIDAQ analysis. Ultimately what I care about is the
%audio in one form or another. Let's just put it here.

%Starting Variables that I need
%An.subject
%An.run
%An.f0AnaFile
%An.curSess
%An.bTf0b

%An.sRate   % 8000Hz (NIDAQ), 16000Hz (Audapter)
%An.audioM  % Expects a matrix (samples x trials)
%An.audioH  % Expects a matrix (samples x trials)

%An.pertIdx   % Vector of indices pertaining to perturbed trials
%An.contIdx   % Vector of indices pertaining to control trials
%An.pertTrig  % Matrix of start/stop indices pertaining to pert period
%An.contTrig  % Matrix of start/stop indices pertaining to 'pert' period

if isempty(varargin)
    iRF    = 0; 
    f0Flag = 0;
elseif length(varargin) == 1
    iRF    = varargin{1};
    f0Flag = 0;
else
    iRF    = varargin{1};
    f0Flag = varargin{2};
end

%Instatiate the variables we intend to use. 
An = initAudVar(An);

if AudFlag == 1
    %Set some frequency analysis variables
    An.fV = setFreqAnalVar(An.sRate);
    
    %Main script that does the Signal Frequency Analysis
    dirs.audiof0AnalysisFile = fullfile(dirs.SavResultsDir, An.f0AnaFile);

    %Sometimes frequency analysis takes a while, this allows you to save
    %results from last time if you want to. 
    if exist(dirs.audiof0AnalysisFile, 'file') == 0 || f0Flag == 1
        ET = tic;
        f0A = signalFrequencyAnalysis(dirs, An.fV, An.audioPP, An.expTrigsPP, An.pertIdx, An.contIdx, An.bTf0b, 2);
%         [f0A.time_audio, f0A.audioHf0, f0A.fsA] = signalFrequencyAnalysis(dirs, An.fV, An.audioH, An.bTf0b, 1);
        
        elapsed_time = toc(ET)/60;
        fprintf('\nElapsed Time: %f (min)\n', elapsed_time)
        save(dirs.audiof0AnalysisFile, 'f0A')
    else
        load(dirs.audiof0AnalysisFile)
    end

    An.fsA     = f0A.fsA;
    An.f0b     = f0A.f0b;
    An.audiof0 = f0A.audiof0;
    An.audiof0Pert = f0A.audiof0Pert;
    An.audiof0Cont = f0A.audiof0Cont;
%     An.audioMf0   = f0A.audioMf0;
%     An.audioHf0   = f0A.audioHf0;

    %Smooth the f0 data
%     An.audiof0S   = smoothf0(An.audiof0);
%     An.audioHf0S   = smoothf0(An.audioHf0);

%     %Normalize f0 and convert to cents
%     prePert       = (0.5 < An.time_audio & 1.0 > An.time_audio);
%     An.trialf0b   = mean(An.audioMf0S(prePert,:),1);
%     An.f0b        = mean(An.trialf0b);
% 
%     An.audioMf0_norm = normf0(An.audioMf0S, An.trialf0b);
%     An.audioHf0_norm = normf0(An.audioHf0S, An.trialf0b);
% 
%     %Find the Perturbed Trials
%     An.audioMf0_p = parseTrialTypes(An.audioMf0_norm, An.pertIdx);
%     An.audioHf0_p = parseTrialTypes(An.audioHf0_norm, An.pertIdx);
%     An.audioMf0_c = parseTrialTypes(An.audioMf0_norm, An.contIdx);
%     An.audioHf0_c = parseTrialTypes(An.audioHf0_norm, An.contIdx);
% 
%     %Find troublesome trials and remove
%     [An.audioMf0_pPP, An.audioHf0_pPP, An.pertTrigPP, An.numPertTrialsPP] = audioPostProcessing(An.time_audio, An.audioMf0_p, An.audioHf0_p, An.pertTrig, An.curSess, 'Pert');
%     [An.audioMf0_cPP, An.audioHf0_cPP, An.contTrigPP, An.numContTrialsPP] = audioPostProcessing(An.time_audio, An.audioMf0_c, An.audioHf0_c, An.contTrig, An.curSess, 'Cont');

    [An.secTime, An.audioMf0_Secp] = sectionAudioData(An.audiof0Pert);

    %Section the data around onset and offset
%     [An.secTime, An.audioMf0_Secp] = sectionAudioData(An.time_audio, An.audioMf0_pPP, An.pertTrigPP);
%     [An.secTime, An.audioHf0_Secp] = sectionAudioData(An.time_audio, An.audioHf0_pPP, An.pertTrigPP);
%     [An.secTime, An.audioMf0_Secc] = sectionAudioData(An.time_audio, An.audioMf0_cPP, An.contTrigPP);
%     [An.secTime, An.audioHf0_Secc] = sectionAudioData(An.time_audio, An.audioHf0_cPP, An.contTrigPP);

    %Mean around the onset and offset
    An.audioMf0_meanp = meanAudioData(An.audioMf0_Secp);
    An.audioHf0_meanp = meanAudioData(An.audioHf0_Secp);
    An.audioMf0_meanc = meanAudioData(An.audioMf0_Secc);
    An.audioHf0_meanc = meanAudioData(An.audioHf0_Secc); 

    %The Inflation Response
    if iRF == 1
        [An.respVar, An.respVarMean, An.respVarSD, An.InflaStimVar] = InflationResponse(An.secTime, An.audioMf0_Secp);
    end
end
end

function An = initAudVar(An)
%Initialize some variables to keep track of them

An.time_audio     = []; %time vector of audio samples recorded
An.fsA            = []; %sampling rate of audio samples
An.audioMf0       = []; %Raw Microphone Audio Data
An.audioHf0       = []; %Raw Headphone Audio Data
An.audioMf0S      = [];
An.audioHf0S      = [];
An.trialf0b       = []; %Per Trial calculated f0
An.f0b            = []; %Average trial f0
An.audioMf0_norm  = []; %Normalized mic data
An.audioHf0_norm  = []; %Normalized head data
An.audioMf0_p     = []; %Perturbed mic data (norm)
An.audioHf0_p     = []; %Pertrubed head data (norm)
An.audioMf0_c     = []; %Control mic data (norm)
An.audioHf0_c     = []; %Control head data (norm)
An.audioMf0_pPP   = []; %
An.audioHf0_pPP   = [];
An.numPertTrialsPP = []; 
An.pertTrigPP      = [];
An.audioMf0_cPP    = []; 
An.audioHf0_cPP    = [];
An.numContTrialsPP = [];
An.contTrigPP      = [];

An.secTime        = [];
An.audioMf0_Secp  = []; 
An.audioHf0_Secp  = [];
An.audioMf0_Secc  = []; 
An.audioHf0_Secc  = [];

An.audioMf0_meanp = []; 
An.audioHf0_meanp = [];
An.audioMf0_meanc = []; 
An.audioHf0_meanc = [];

An.respVar        = []; 
An.respVarMean    = [];
An.respVarSD      = [];
An.InflaStimVar   = [];
end

function fV = setFreqAnalVar(sRate)

%Identify a few analysis varaibles
fV.sRate      = sRate;

fV.freqCutOff = 400;
fV.win        = 0.015;      % seconds
fV.fsA        = 1/fV.win;
fV.winP       = fV.win*sRate;
fV.pOV        = 0.80;       % 80% overlap
fV.tStepP     = fV.winP*(1-fV.pOV);
end

function [f0A] = signalFrequencyAnalysis(dirs, fV, audioPP, expTrigs, pertIdx, contIdx, bTf0b, flag)
[~, numTrial] = size(audioPP);

fs = fV.sRate;

audiof0 = []; allf0b = [];
for j = 1:numTrial %Trial by Trial 
    time     = audioPP(j).time;
    mic      = audioPP(j).mic;
    head     = audioPP(j).head;
    expTrig  = expTrigs(j,:);

    numSamp = length(mic);

    trialWinSt   = round(1:fV.tStepP:(numSamp-fV.winP));  % Window start frames based on length of voice onset mic
    numWinStT    = length(trialWinSt);                    % Number of windows based on WinSt

    timef0 = []; micf0 = []; headf0 = [];
    for i = 1:numWinStT
        winIdx  = trialWinSt(i):trialWinSt(i)+ fV.winP - 1;

        timeWin = mean(time(winIdx));            
        micWin  = mic(winIdx);
        headWin = head(winIdx);

        % What is the f0??
%             f0Win = dfCalcf0Chile(voiceWinHP, fs);
%             if isnan(f0Win); f0Win = bTf0b; end

        micf0Win  = simpleAutoCorr(micWin, fV);
        headf0Win = simpleAutoCorr(headWin, fV);

        timef0 = cat(1, timef0, timeWin);
        micf0  = cat(1, micf0, micf0Win);
        headf0 = cat(1, headf0, headf0Win);
    end

    %Smooth the trace over 10 samples
    micf0S  = smooth(micf0, 10);
    headf0S = smooth(headf0, 10);

    onSetInd   = find(timef0 <= expTrig(1));
    onsetI     = onSetInd(end);
    prePertPer = 1:onsetI;  % Everything is already cut starting at 0.5s before pert onset
    trialf0b   = mean(micf0S(prePertPer));

    micf0Norm  = normf0(micf0S, trialf0b);
    headf0Norm = normf0(headf0S, trialf0b);

    gordon = find(micf0Norm > 500 | micf0Norm < -500);
    if ~isempty(gordon)
        disp('Mamma Mia')
    end

    audiof0(j).timef0 = timef0;
    audiof0(j).micf0  = micf0S;  
    audiof0(j).headf0 = headf0S;
    audiof0(j).trialf0b   = trialf0b;
    audiof0(j).micf0Norm  = micf0Norm;
    audiof0(j).headf0Norm = headf0Norm; 
    audiof0(j).expTrig    = expTrig;
    allf0b = cat(1, allf0b, trialf0b);
end
fsA = fV.fsA;

audiof0Pert = audiof0(pertIdx);
audiof0Cont = audiof0(contIdx);

f0b = mean(allf0b);

f0A.fsA         = fsA;
f0A.f0b         = f0b;
f0A.audiof0     = audiof0;
f0A.audiof0Pert = audiof0Pert;
f0A.audiof0Cont = audiof0Cont;
end

function f0Win = simpleAutoCorr(voice, fV)
%Simple version of an autocorrelation for finding pitch

fs            = fV.sRate;
win           = fV.winP;
[autoC, lags] = autocorr(voice, win-1);
[pks, pkInd]  = findpeaks(autoC);
if isempty(pks)
    FLag      = lags(end);
else
    [~, mxInd]= max(pks);
    FLag      = lags(pkInd(mxInd));
end
per           = FLag/fs;
f0Win         = 1/per;
end

function audioS = smoothf0(audio)
[~, numTrial] = size(audio);

audioS = [];
for ii = 1:numTrial
    audioSmooth = smooth(audio(:,ii), 10);
    audioS      = cat(2, audioS, audioSmooth);
end
end

function audio_norm = normf0(audio, f0b)
%audio_norm = normf0(audio, f0b) is a function that takes a set of audio signals and 
%normalizes each piece of audio by its corresponding fundamental fequency.
%It is expected that audio and f0b will be a matrices of size numSamp x
%numTrial 

[~, numTrial] = size(audio);

audio_norm = [];
for ii = 1:numTrial
    audio_trial = 1200*log2(audio(:,ii)./f0b(ii));
    audio_norm  = cat(2, audio_norm, audio_trial);
end
end

function signalParse = parseTrialTypes(signal, idx)
%Expects trials to be in columns 

signalParse = signal(:, idx); %This is a little lazy I know. Get over it. 
end

function [audioNormMPP, audioNormHPP, trigsPP, numTrialTypePP] = audioPostProcessing(time, audioNormM, audioNormH, trigs, curSess, type)
%This function checks to see if there are any realllllly weird f0 values as
%a result of the spectral analysis. This throws away trials and tells you
%when it happens. It combines new data audio files of the saved trials, and
%ammends the trig matrix of those files removed. Also gives a new value of
%total number of trials kept.
[~, numTrialType] = size(audioNormM);

timeInd = find(time > 0.5 & time < 3.5);

audioNormMPP   = [];
audioNormHPP   = [];
trigsPP        = [];
numTrialTypePP = 0; 
for ii = 1:numTrialType
    ind = find(audioNormM(timeInd,ii) >= 500 | audioNormM(timeInd,ii) <=  -500);
    if ~isempty(ind)
        fprintf('Threw away %s %s trial %s\n', curSess, type, num2str(ii))
    else
        trigsPP        = cat(1, trigsPP, trigs(ii,:));
        audioNormMPP   = cat(2, audioNormMPP, audioNormM(:,ii));
        audioNormHPP   = cat(2, audioNormHPP, audioNormH(:,ii));
        numTrialTypePP = numTrialTypePP + 1;
    end
end
end

function [secTime, secAudio] = sectionAudioData(audiof0)
[~, numTrial] = size(audiof0);
preEve  = 0.5; posEve = 1.0;
% per     = 1/fs;
% preEveP = preEve*fs;
% posEveP = posEve*fs-1;
% trigsR   = round2matchfs(trigs);

secAudio   = [];
OnsetSecs  = [];
OffsetSecs = [];
for ii = 1:numTrial
    time = audiof0(ii).timef0;
    micf0 = audiof0(ii).micf0Norm;
    headf0 = audiof0(ii).headf0Norm;
    trigs  = audiof0(ii).expTrig;
    
    OnsetT   = trigs(1);
    OffsetT  = trigs(2);
    
    OnsetTSt = round(OnsetT - preEve, 4);   % Accurate to nearest ms
    OnsetTSp = round(OnsetT + posEve, 4);   % Accurate to nearest ms
    OnsetSpan = time >= OnsetTSt & time <= OnsetTSp;
    
    OffsetTSt = round(OffsetT - preEve, 4); % Accurate to nearest ms
    OffsetTSp = round(OffsetT + posEve, 4); % Accurate to nearest ms
    OffsetSpan = time >= OffsetTSt & time <= OffsetTSp;
    
%     OnsetI  = find(time == trigsR(ii,1));
%     OffsetI = find(time == trigsR(ii,2));
%     
%     OnsetISt = OnsetI - preEveP;
%     OnsetISp = OnsetI + posEveP;
%     
%     OffsetISt = OffsetI - preEveP;
%     OffsetISp = OffsetI + posEveP;
        
    OnsetSec  = micf0(OnsetSpan);
    OffsetSec = micf0(OffsetSpan);
    
    OnsetSecs  = cat(2, OnsetSecs, OnsetSec);
    OffsetSecs = cat(2, OffsetSecs, OffsetSec);
end

numSamp = length(OnsetSec);

secTime = linspace(-preEve, posEve, numSamp);
secAudio(:,:,1) = OnsetSecs; 
secAudio(:,:,2) = OffsetSecs;
end

function y = round2matchfs(x)
%This expects a decimal number as input
%Input can be given as a set

y = round(x.*200)./200;
end

function meanAudio = meanAudioData(secAudio)

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

function [respVar, respVarm, respVarSD, InflaStimVar] = InflationResponse(secTime, secAudio)
[L, numTrial, ~] = size(secAudio); %Only look at Onsets

ir.numTrial = numTrial;
ir.time     = secTime;
ir.iAtOnset = find(secTime == 0); %Ind
ir.tAtOnset = 0;
ir.vAtOnset = [];
ir.iPostOnsetR = find(0 <= secTime & .20 >= secTime); %Ind
ir.iAtMin = [];
ir.tAtMin = [];
ir.vAtMin = [];
ir.stimMag = [];
ir.iAtResp = L; %the last ind
ir.tAtResp = ir.time(L);
ir.vAtResp = [];
ir.respMag = [];
ir.respPer = [];

shpInds = [];
tAtMin  = []; stimMag = [];
respMag = []; respPer = [];
for i = 1:numTrial
    onset = secAudio(:,i,1); %First depth dim in Onset
    ir.vAtOnset = onset(ir.iAtOnset);

    [minOn, minIdx] = min(onset(ir.iPostOnsetR));
    ir.iAtMin = ir.iPostOnsetR(minIdx);
    ir.tAtMin = ir.time(ir.iAtMin);
    ir.vAtMin = minOn;
    ir.stimMag = ir.vAtMin - ir.vAtOnset;
    
    ir.vAtResp = onset(ir.iAtResp);
    ir.respMag = ir.vAtResp - ir.vAtMin;
    
    ir.respPer = 100*(ir.respMag/abs(ir.stimMag));
    
    if ir.stimMag == 0
        ir.respPer = 0.0;
    end
    
%     subplot(2,5,i)
%     plot(secTime, onset)
    
%     shpInd   = [ir.iAtOnset ir.iAtMin ir.iAtResp];  
%     shpInds  = cat(1, shpInds, shpInd); 
    
    tAtMin   = cat(1, tAtMin, ir.tAtMin);
    stimMag  = cat(1, stimMag, ir.stimMag); 
    respMag  = cat(1, respMag, ir.respMag); 
    respPer  = cat(1, respPer, ir.respPer);
end

respVar  = [tAtMin stimMag respMag respPer];
respVarm = mean(respVar, 1);
respVarSD = std(respVar, 0, 1);

InflaStimVar = [respVar(1) respVarm(2)];
end