function audiof0 = dfAnalysisAudioQuick(DRF)

close all
rawData  = DRF.rawData;
expParam = DRF.expParam;
expParam.frameLenDown    = expParam.frameLen/expParam.downFact;

% Find the indices at which voicing starts
[voiceInd] = preProcessVoice(rawData, expParam.frameLenDown);

fV = setFreqAnalVar(expParam.sRateAnal, voiceInd);

% Some quick pitch analysis of each trial. 
[audiof0] = signalFrequencyAnalysis(fV, rawData);

% plotBaseTrials(audiof0)
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

function [audiof0] = signalFrequencyAnalysis(fV, rawData)
[numTrial, ~] = size(rawData);

fs   = fV.sRate;

%Low-Pass filter for the given cut off frequency
[B,A]    = butter(4,(fV.freqCutOff)/(fs/2));

audiof0 = [];
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
         
        f0Win = quickDirtyAutoCorr(voiceWinHP, fV);      
        
        timef0  = cat(1, timef0, timeWin);
        voicef0 = cat(1, voicef0, f0Win);
    end
    
    audiof0(j).time = smooth(timef0,10);
    audiof0(j).f0   = smooth(voicef0,10);  
end
end

function f0Win = quickDirtyAutoCorr(voice, fV)
fs            = fV.sRate;
[autoC, lags] = autocorr(voice, fV.winP-1);
[pks, pkInd]  = findpeaks(autoC);
[~, mxInd]    = max(pks);
FLag          = lags(pkInd(mxInd));
per           = FLag/fs;
f0Win         = 1/per;
end

function plotBaseTrials(audiof0)
numTrial = length(audiof0);

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
    title(['Trial ' num2str(j)]);
    axis([0 4 200 240])
    box off
    
    set(gca,'FontSize', 11,...
        'FontWeight','bold')
    
end

suptitle('Pilot28')
end