function res = dfAnalysisAudioQuick(DRF, varargin)

if isempty(varargin)
    pltFlg = 0;
else
    pltFlg = varargin{1};    
end

rawData  = DRF.rawData;
expParam = DRF.expParam;
subj     = expParam.subject;
run      = expParam.run;
rmsB     = expParam.rmsB;
frameLenDown    = expParam.frameLen/expParam.downFact;

% Find the indices at which voicing starts
[voiceInd] = preProcessVoice(rawData, frameLenDown);

fV = setFreqAnalVar(expParam.sRateAnal, voiceInd);

% Some quick pitch analysis of each trial. 
[audiof0, trialf0, f0Bounds, audioRMS, elapsed_time] = signalFrequencyAnalysis(fV, rawData, rmsB);

res.audiof0  = audiof0;
res.trialf0  = trialf0;
res.meanf0   = mean(trialf0);
res.f0Bounds = f0Bounds;
res.audioRMS = audioRMS;
res.meanRMS  = mean(audioRMS);
res.anaTime  = round(elapsed_time, 2);

if pltFlg == 1
    plotBaseTrials(subj, run, res)
end
end

function [voiceInd] = preProcessVoice(rawData, frameLen)
[numTrial, ~] = size(rawData);

thresh = 0.50;
voiceInd = [];
for jj = 1:numTrial
    trialData = rawData(jj);
    trialDataRMS = trialData.rms(:,1);
    maxtrialRMS  = max(trialDataRMS);
    ind = find(trialDataRMS > maxtrialRMS*thresh);
    indSig = ind(1)*frameLen;
    voiceInd = cat(1, voiceInd, indSig);    
end
end

function fV = setFreqAnalVar(sRate, voiceInd)

%Identify a few analysis varaibles
fV.sRate      = sRate;

fV.freqCutOff = 400;
fV.win        = 0.050;       % seconds
fV.fsA        = 1/fV.win;
fV.winP       = fV.win*sRate;
fV.pOV        = 0.50;        % 80% overlap
fV.tStepP     = fV.winP*(1-fV.pOV);
fV.voiceInd   = voiceInd;
end

function [audiof0, trialf0, f0Bounds, audioRMS, elapsed_time] = signalFrequencyAnalysis(fV, rawData, rmsB)
ET = tic;
[numTrial, ~] = size(rawData);

fs   = fV.sRate;

audiof0 = []; trialf0 = []; audioRMS = [];
f0Max = 0; f0Min = 0;
for j = 1:numTrial %Trial by Trial
    trialData = rawData(j);     % Grab the data set for this trial
    voiceInd  = fV.voiceInd(j); % Pre-calculated RMS threshold
   
    voiceW    = trialData.signalIn;   % Microphone Channel
    numSampW  = length(voiceW);       % Length of Mic Channel
    voiceT    = voiceW(voiceInd:end); % Microphone Channel, starting at Voice Onset
    numSampT  = length(voiceT);       % Length of Mic Channel, starting at Voice Onset
    
    timeW     = (0:1/fs:(numSampW-1)/fs)';           % Time vector for full mic
    timeT     = timeW(voiceInd:voiceInd+numSampT-1); % Time vector for mic starting at Voice Onset
    
    trialWinSt = round(1:fV.tStepP:(numSampT-fV.winP)); % Window start frames based on length of voice onset mic
    numWinStT  = length(trialWinSt);                    % Number of windows based on WinSt
    
    timef0 = []; voicef0 = [];
    for i = 1:numWinStT
%         fprintf('Trial %d, window %d\n', j, i)
        winIdx  = trialWinSt(i):trialWinSt(i)+ fV.winP - 1;
        timeWin = mean(timeT(winIdx));
        
        % What is the f0??
        voiceWin   = voiceT(winIdx);
        f0Win = simpleAutoCorr(voiceWin, fV);      
        
        timef0  = cat(1, timef0, timeWin);
        voicef0 = cat(1, voicef0, f0Win);
    end
    
    voicef0S = smooth(voicef0, 10);
    
    audiof0(j).time = timef0;
    audiof0(j).f0   = voicef0S;
    
    %Overall steady-state f0
    voicef0Mean = round(mean(voicef0S), 2);
    trialf0 = cat(1, trialf0, voicef0Mean);
    
    %Range of frequencies in this trial
    trialf0Max = max(voicef0S(10:end));
    trialf0Min = min(voicef0S(10:end));
    if j == 1 || trialf0Max > f0Max
        f0Max = trialf0Max;
    end

    if j == 1 || trialf0Min < f0Min
        f0Min = trialf0Min;
    end
    
    %RMS Calculation
    voiceRMS = dfCalcMeanRMS(trialData, rmsB);
    audioRMS = cat(1, audioRMS, voiceRMS);    
end

f0Max = round(f0Max) + 5;
f0Min = round(f0Min) - 5;
f0Bounds = [f0Min f0Max];

elapsed_time = toc(ET)/60;
fprintf('\nElapsed Time: %f (min)\n', elapsed_time)
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

function plotBaseTrials(subj, run, res)
audiof0 = res.audiof0;
trialf0 = res.trialf0;
meanf0  = res.meanf0;
lB      = res.f0Bounds(1);
hB      = res.f0Bounds(2);
anaTime = res.anaTime;

numTrial = length(audiof0);
if numTrial > 3; numTrial = 3; end

plotpos = [10 30];
plotdim = [600 900];
pitchTrace = figure('Color', [1 1 1]);
set(pitchTrace, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

ha = tight_subplot(numTrial, 1, [0.1 0.05],[0.08 0.05],[0.1 0.05]);

for j = 1:numTrial
    axes(ha(j))
    plot(audiof0(j).time, audiof0(j).f0, 'b', 'LineWidth', 2)
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title(['Trial ' num2str(j) ' f0: ' num2str(trialf0(j)) ' Hz']);
    axis([0 4 lB hB])
    box off
    
    set(gca,'FontSize', 11,...
            'FontWeight','bold')
end
suptitle({[subj run '  f0Ave: ' num2str(meanf0) ' Hz'],...
          ['Analysis Time: ' num2str(anaTime) ' sec']});
end