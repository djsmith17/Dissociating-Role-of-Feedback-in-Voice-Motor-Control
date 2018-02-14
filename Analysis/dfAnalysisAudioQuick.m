function quickResult = dfAnalysisAudioQuick(DRF, varargin)

if isempty(varargin)
    pltFlg = 0;
else
    pltFlg = varargin{1};    
end

rawData  = DRF.rawData;
expParam = DRF.expParam;
subj     = expParam.subject;
expParam.frameLenDown    = expParam.frameLen/expParam.downFact;

% Find the indices at which voicing starts
[voiceInd] = preProcessVoice(rawData, expParam.frameLenDown);

fV = setFreqAnalVar(expParam.sRateAnal, voiceInd);

% Some quick pitch analysis of each trial. 
[audiof0, trialf0, f0Bounds, audioRMS] = signalFrequencyAnalysis(fV, rawData);

quickResult.audiof0  = audiof0;
quickResult.trialf0  = trialf0;
quickResult.meanf0   = mean(trialf0);
quickResult.f0Bounds = f0Bounds;
quickResult.audioRMS = audioRMS;
quickResult.meanRMS  = mean(audioRMS);

if pltFlg == 1
    plotBaseTrials(subj, quickResult)
end
end

function [voiceInd] = preProcessVoice(rawData, frameLen)
[numTrial, ~] = size(rawData);

thresh = 0.01;
voiceInd = [];
for jj = 1:numTrial
    trialData = rawData(jj);
    ind = find(trialData.rms(:,1) > thresh);
    indSig = ind(1)*frameLen;
    voiceInd = cat(1, voiceInd, indSig);    
end
end

function fV = setFreqAnalVar(sRate, voiceInd)

%Identify a few analysis varaibles
fV.sRate      = sRate;

fV.freqCutOff = 400;
fV.win        = 0.015;       % seconds
fV.fsA        = 1/fV.win;
fV.winP       = fV.win*sRate;
fV.pOV        = 0.80;        % 80% overlap
fV.tStepP     = fV.winP*(1-fV.pOV);
fV.voiceInd   = voiceInd;
end

function [audiof0, trialf0, f0Bounds, audioRMS] = signalFrequencyAnalysis(fV, rawData)
[numTrial, ~] = size(rawData);

fs   = fV.sRate;

%Low-Pass filter for the given cut off frequency
[B,A]    = butter(4,(fV.freqCutOff)/(fs/2));

audiof0 = []; trialf0 = []; audioRMS = [];
f0Max = 0; f0Min = 0;
for j = 1:numTrial %Trial by Trial
    trialData = rawData(j);
    voiceInd  = fV.voiceInd(j);
   
    voiceW    = trialData.signalIn;
    numSampW  = length(voiceW);
    voiceT    = voiceW(voiceInd:end);
    numSampT  = length(voiceT);
    
    timeW     = (0:1/fs:(numSampW-1)/fs)';
    timeT     = timeW(voiceInd:voiceInd+numSampT-1);       
    
    trialWinSt = round(1:fV.tStepP:(numSampT-fV.winP));
    numWinStT  = length(trialWinSt);
    
    timef0 = []; voicef0 = [];
    for i = 1:numWinStT
        winIdx  = trialWinSt(i):trialWinSt(i)+ fV.winP - 1;
        timeWin = mean(timeT(winIdx));
        
        voiceWin   = voiceT(winIdx);
        voiceWinHP = filtfilt(B, A, voiceWin);
         
        f0Win = simpleAutoCorr(voiceWinHP, fV);      
        
        timef0  = cat(1, timef0, timeWin);
        voicef0 = cat(1, voicef0, f0Win);
    end
    
    audiof0(j).time = smooth(timef0,10);
    audiof0(j).f0   = smooth(voicef0,10);
    
    %Overall steady-state f0
    voicef0Mean = round(mean(voicef0), 2);
    trialf0 = cat(1, trialf0, voicef0Mean);
    
    %Range of frequencies in this trial
    trialf0Max = max(voicef0(10:end));
    trialf0Min = min(voicef0(10:end));
    if j == 1 || trialf0Max > f0Max
        f0Max = trialf0Max;
    end

    if j == 1 || trialf0Min < f0Min
        f0Min = trialf0Min;
    end
    
    %RMS Calculation
    voiceRMS = dfCalcMeanRMS(trialData);
    audioRMS = cat(1, audioRMS, voiceRMS);    
end

f0Max = round(f0Max) + 5;
f0Min = round(f0Min) - 5;
f0Bounds = [f0Min f0Max];
end

function f0Win = simpleAutoCorr(voice, fV)
%Simple version of an autocorrelation for finding pitch

fs            = fV.sRate;
win           = fV.winP;
[autoC, lags] = autocorr(voice, win-1);
[pks, pkInd]  = findpeaks(autoC);
[~, mxInd]    = max(pks);
FLag          = lags(pkInd(mxInd));
per           = FLag/fs;
f0Win         = 1/per;
end

function plotBaseTrials(subj, quickResult)
audiof0 = quickResult.audiof0;
trialf0 = quickResult.trialf0;
meanf0 = quickResult.meanf0;
lB = quickResult.f0Bounds(1);
hB = quickResult.f0Bounds(2);

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
suptitle([subj '  f0Ave: ' num2str(meanf0) ' Hz']);
end