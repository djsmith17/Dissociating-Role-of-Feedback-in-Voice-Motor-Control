function dfAnalysisAudioQuick(DRF)

rawData  = DRF.rawData;
expParam = DRF.expParam;
expParam.numSamp            = expParam.sRateAnal*expParam.trialLen;
expParam.frameLenDown       = expParam.frameLen/expParam.downFact;

% Find the indices at which voicing starts
[voiceInd] = preProcessVoice(rawData, expParam.frameLenDown);

fV = setFreqAnalVar(expParam.sRateAnal, expParam.numSamp, voiceInd);

% Some quick pitch analysis of each trial. 
[audiof0] = signalFrequencyAnalysis(fV, rawData);

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

function fV = setFreqAnalVar(sRate, numSamp, voiceInd)

%Identify a few analysis varaibles
fV.sRate      = sRate;
fV.numSamp    = numSamp;
fV.time       = (0:1/sRate:(numSamp-1)/sRate)';

fV.freqCutOff = 500;
fV.win        = 0.05;       % seconds
fV.fsA        = 1/fV.win;
fV.winP       = fV.win*sRate;
fV.pOV        = 0.60;        % 60% overlap
fV.tStepP     = fV.winP*(1-fV.pOV);
fV.voiceInd   = voiceInd;
fV.winSts     = 1:fV.tStepP:(numSamp-fV.winP);
fV.numWin     = length(fV.winSts);
end

function [audiof0] = signalFrequencyAnalysis(fV, rawData)
[numTrial, ~] = size(rawData);

time = fV.time;
fs   = fV.sRate;

%Low-Pass filter for the given cut off frequency
[B,A]    = butter(4,(fV.freqCutOff)/(fs/2));

audiof0 = [];
for j = 1:numTrial %Trial by Trial
    trialData   = rawData(j);
    trialVoiceW = trialData.signalIn(1:fV.numSamp);
    
    voiceInd     = fV.voiceInd(j);
    trialTime    = time(voiceInd:end);
    trialVoice   = trialVoiceW(voiceInd:end);
    numSampTrial = length(trialVoice);    
    
    trialWinSt = 1:fV.tStepP:(numSampTrial-fV.winP);
    numWinStT  = length(trialWinSt);
    
    timef0 = []; voicef0 = [];
    for i = 1:numWinStT
        winIdx  = trialWinSt(i):trialWinSt(i)+ fV.winP - 1;
        timeWin = mean(trialTime(winIdx));
        
        voiceWin   = trialVoice(winIdx);
        voiceWinHP = filtfilt(B, A, voiceWin);
         
%         f0Win = quickDirtyAutoCorr(voiceWinHP, fV)
        f0Win = quickFFT(voiceWin, fs, fV);
%         f0Win = dfCalcf0Chile(voiceWinHP, fs);       
        
        timef0  = cat(1, timef0, timeWin);
        voicef0 = cat(1, voicef0, f0Win);
    end
    
    audiof0(j).time = timef0;
    audiof0(j).f0   = voicef0;    
end
end

function f0Win = quickDirtyAutoCorr(voice, fV)
[autoC, lags] = autocorr(voice, fV.winP-1);
[~, pkInd]    = findpeaks(autoC);
FLag          = lags(pkInd(1));
per           = FLag/fs;
f0Win         = 1/per;
end

function f0 = quickFFT(voice, fs, fV)

L = length(voice);
NFFT = 2^nextpow2(L);

winN = fV.winP;
nOverLap = winN*fV.pOV;

[pxx, f] = pwelch(voice, winN, nOverLap, NFFT, fs);

[~, ind] = max(pxx);
f0 = f(ind);
end