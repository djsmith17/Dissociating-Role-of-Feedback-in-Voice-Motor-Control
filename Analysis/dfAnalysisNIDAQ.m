function [niAn, res] = dfAnalysisNIDAQ(dirs, expParam, DAQin)
%A quick reference
%
%Pert: Perturbation signal
%P:    Pressure sensor signal
%FC:   Force Sensor Collar
%FN:   Force Sensor Neck
%
%Trig: Trigger values where onset and offset occur
%DN:   Down Sampled (and smoothed)

fprintf('\nStarting NIDAQ Analysis\n')

sv2File = 0;
[r, c, n] = size(DAQin);
sRate = expParam.sRateQ;

%Identify some starting variables
niAn.subject  = expParam.subject;
niAn.run      = expParam.run;
if isfield(expParam, 'curSess')
    niAn.curSess  = expParam.curSess;  %Short hand of experiment details
else
    niAn.curSess  = [expParam.subject expParam.run]; 
end
niAn.trialType    = expParam.trialType;

niAn.dnSamp   = 10;
niAn.sRate    = sRate;
niAn.numSamp  = r;
niAn.numTrial = n;
niAn.numCh    = c;
niAn.expTrigs = expParam.trigs(:,:,1); %Time

%Find all the pertrubed trials
[niAn.PertTrials, niAn.pertIdx] = find(niAn.trialType == 1);
niAn.numPertTrials = sum(niAn.PertTrials);

%Identify a few analysis varaibles
niAn.win    = 0.05;  %seconds
niAn.winP   = niAn.win*niAn.sRate;
niAn.pOV    = 0.60;  %60% overlap
niAn.tStepP = niAn.winP*(1-niAn.pOV);
niAn.winSts = 1:niAn.tStepP:(niAn.numSamp-niAn.winP);
niAn.numWin = length(niAn.winSts);
niAn.freqCutOff = 400;

%Unpack the NIDAQ raw data set
niAn.sRateDN  = sRate/niAn.dnSamp;
niAn.time     = (0:1/sRate:(r-1)/sRate)';
niAn.pertSig  = squeeze(DAQin(:,1,:));
niAn.sensorFC = squeeze(DAQin(:,2,:));
niAn.sensorFN = squeeze(DAQin(:,3,:));
niAn.sensorP  = squeeze(DAQin(:,4,:));
niAn.audioM   = squeeze(DAQin(:,5,:));
niAn.audioH   = squeeze(DAQin(:,6,:));
niAn.sensorO  = squeeze(DAQin(:,7,:));

%Preprocessing and downsampling
niAn.sensorFC_aug = 4*(niAn.sensorFC-2);
niAn.sensorFN_aug = 4*(niAn.sensorFN-2);

[B,A] = butter(4, 10/(sRate/2)); %Low-pass filter under 10
niAn.sensorFC = filter(B,A,abs(niAn.sensorFC));
niAn.sensorFN = filter(B,A,abs(niAn.sensorFN));

niAn.time_DN     = dnSampleSignal(niAn.time, niAn.dnSamp);    % DownSampled Time
niAn.pertSig_DN  = dnSampleSignal(niAn.pertSig, niAn.dnSamp); % DownSampled Perturbation Signal
niAn.sensorP_DN  = dnSampleSignal(niAn.sensorP, niAn.dnSamp);
niAn.sensorFC_DN = dnSampleSignal(niAn.sensorFC, niAn.dnSamp);
niAn.sensorFN_DN = dnSampleSignal(niAn.sensorFN, niAn.dnSamp);

%Parse out the pertrubed trials
niAn.pertSig_p  = niAn.pertSig_DN(:, niAn.pertIdx);   %Grab all the Catch Trials
niAn.sensorP_p  = niAn.sensorP_DN(:, niAn.pertIdx);   %Grab all the Catch Trials
niAn.sensorFC_p = niAn.sensorFC_DN(:, niAn.pertIdx);  %Grab all the Catch Trials
niAn.sensorFN_p = niAn.sensorFN_DN(:, niAn.pertIdx);  %Grab all the Catch Trials

niAn.expTrigs_p = niAn.expTrigs(niAn.trialType == true, :);
[niAn.pertTrig, niAn.idxPert] = findPertTrigs(niAn.time_DN, niAn.pertSig_p, niAn.sRateDN);
[niAn.presTrig, niAn.idxPres] = findPertTrigs(niAn.time_DN, niAn.sensorP_p, niAn.sRateDN);
[niAn.fSCTrig, niAn.idxFC]    = findPertTrigs(niAn.time_DN, niAn.sensorFC_p, niAn.sRateDN);  
[niAn.fSNTrig, niAn.idxFN]    = findPertTrigs(niAn.time_DN, niAn.sensorFN_p, niAn.sRateDN); 

[niAn.lagsPres, niAn.meanLagTimeP] = calcMeanLags(niAn.pertTrig, niAn.presTrig);
[niAn.lagsFC, niAn.meanLagTimeFC]  = calcMeanLags(niAn.pertTrig, niAn.fSCTrig);
[niAn.lagsFN, niAn.meanLagTimeFN]  = calcMeanLags(niAn.pertTrig, niAn.fSNTrig);

niAn.indPressures   = [];
niAn.rangePressures = [];
niAn.timePressures  = [];
for ii = 1:niAn.numPertTrials
    [endRiseInd, startFallInd] = findCrossings(niAn.sensorP_p(:,ii), niAn.sRateDN);
%     [endRiseInd, startFallInd] = findCrossingsManual(niAn.sensorP_DN(:,ii), niAn.sRateDN);
    
%     [maxP, maxInd] = max(niAn.sensorP_DN(:,ii));
%     [minP  ] = niAn.sensorP_DN(niAn.idxPert(ii,2), ii);

    onsetPressure  = round(100*niAn.sensorP_p(endRiseInd, ii))/100;
    offsetPressure = round(100*niAn.sensorP_p(startFallInd, ii))/100;
    onsetTime  = round(100*niAn.time_DN(endRiseInd))/100;
    offsetTime = round(100*niAn.time_DN(startFallInd))/100;
    
    niAn.indPressures   = cat(1, niAn.indPressures, [endRiseInd, startFallInd]);
    niAn.timePressures  = cat(1, niAn.timePressures, [onsetTime offsetTime]);
    niAn.rangePressures = cat(1, niAn.rangePressures, [onsetPressure offsetPressure]);
end
niAn.meanRangePressure = mean(niAn.rangePressures, 1);

niAn.riseTimeP = niAn.timePressures(:,1) - niAn.presTrig(:,1);
niAn.meanRiseTimeP = mean(niAn.riseTimeP);

niAn.sensorP_Al = alignSensorData(niAn.sRateDN, niAn.idxPert, niAn.sensorP_p);
niAn.time_Al    = (0:1/niAn.sRateDN :(length(niAn.sensorP_Al)-1)/niAn.sRateDN)';

[niAn.time_audio, niAn.audioMf0] = dfCalcf0Praat(dirs, niAn.audioM, niAn.sRate);
[niAn.time_audio, niAn.audioHf0] = dfCalcf0Praat(dirs, niAn.audioH, niAn.sRate);
% niAn.time_audio = dnSampleSmoothSignal(niAn.time, niAn.winP, niAn.numWin, niAn.winSts);
% niAn.audioMf0   = signalFrequencyAnalysis(niAn.audioM, niAn.sRate, niAn.freqCutOff, niAn.numTrial, niAn.numWin, niAn.winSts, niAn.winP);
% niAn.audioHf0   = signalFrequencyAnalysis(niAn.audioH, niAn.sRate, niAn.freqCutOff, niAn.numTrial, niAn.numWin, niAn.winSts, niAn.winP);
prePert         = (0.5 < niAn.time_audio & 1.0 > niAn.time_audio);
niAn.trialf0b   = mean(niAn.audioMf0(prePert,:),1);
niAn.f0b        = mean(niAn.trialf0b);

niAn.audioMf0_norm = normalizef0(niAn.audioMf0, niAn.trialf0b);
niAn.audioHf0_norm = normalizef0(niAn.audioMf0, niAn.trialf0b);
niAn.audioMf0_p = niAn.audioMf0_norm(:, niAn.pertIdx);
niAn.audioHf0_p = niAn.audioHf0_norm(:, niAn.pertIdx);

lims = identifyLimits(niAn);
res  = packResults(niAn, lims);
end

function sensorDN = dnSampleSignal(sensor, dnSamp)
[numSamp, numTrial] = size(sensor);
numSampDN = numSamp/dnSamp;

sensorDN = zeros(numSampDN, numTrial);
for i = 1:numSampDN
    sensorDN(i,:) = mean(sensor((1:dnSamp) + dnSamp*(i-1),:));
end
end

function sensorDN = dnSampleSmoothSignal(sensor, winP, numWin, winSts)

sensorDN = [];
for iSt = 1:numWin
    winIdx = winSts(iSt):winSts(iSt) + winP - 1;
    sensorDN = cat(1, sensorDN, mean(sensor(winIdx, :)));
end
end

function sensorf0 = signalFrequencyAnalysis(sensor, fs, freqCutOff, numTrial, numWin, winSts, winP)

%Low-Pass filter for the given cut off frequency
[B,A]    = butter(4,(freqCutOff)/(fs/2));

sensorf0 = zeros(numWin, numTrial);
for j = 1:numTrial %Trial by Trial
    sensorHP = filtfilt(B,A,sensor(:,j));
    for i = 1:numWin
        winIdx = winSts(i):winSts(i)+ winP - 1;
        sensorf0(i,j) = dfCalcf0Chile(sensorHP(winIdx), fs);
    end
end
end

function audio_norm = normalizef0(audio, f0b)
[~, numTrial] = size(audio);

audio_norm = [];
for ii = 1:numTrial
    audio_trial = 1200*log2(audio(:,ii)./f0b(ii));
    audio_norm  = cat(2, audio_norm, audio_trial);
end
end

function [trigs, idx] = findPertTrigs(time, sensor, fs)
[numSamp, numTrial] = size(sensor);

baselineIdx = 1:fs*1;

trigs = [];
idx   = [];
for i = 1:numTrial
    m = [];
    for j = 2:numSamp
        m(j) = sensor(j,i) - sensor((j-1),i);
    end
    m(baselineIdx) = 0;
    threshUp = 0.5*max(m);
    threshDn = 0.5*min(m);
    ups = find(m > threshUp);
    dns = find(m < threshDn);
    
    idxSt = ups(1); 
    idxSp = dns(1);       
    trigSt = round(1000*time(idxSt))/1000;
    trigSp = round(1000*time(idxSp))/1000;

    trigs = cat(1, trigs, [trigSt trigSp]);
    idx   = cat(1, idx, [idxSt idxSp]);
end
end

function [lags, lagMeans] = calcMeanLags(pertTrig, sensorTrig)

lags = sensorTrig - pertTrig;
lagsMean = mean(lags, 1);
lagsSTD  = std(lags, 0, 1);

SEM = lagsSTD/sqrt(length(lags)); 
CIM = 1.96*SEM;

lagMeans = [lagsMean, CIM];
end

function sensorAl = alignSensorData(sRate, idx, sensor)
numTrial = length(idx);

sensorAl = [];
for ii = 1:numTrial
    St = idx(ii,1) - sRate*1;
    Sp = idx(ii,1) + sRate*2.5;
    
    sensorAl = cat(2, sensorAl, sensor(St:Sp, ii));
end
end

function [endRiseInd, startFallInd] = findCrossings(sensor, fs)

[B, A] = butter(8, (50/(fs/2)), 'low'); 
sensorFilt = filtfilt(B,A, sensor);

sensDiff = [0; diff(sensorFilt)]*20;
sensDiff2= [0; diff(sensDiff)]*20;

incInds = find(sensDiff > 0.05);
decInds = find(sensDiff < -0.05);

incNonIncre = find((diff(incInds) == 1) == 0); %The indices where the values stop increasing
endRiseInd = incInds(incNonIncre(end));

decNonIncre = find((diff(decInds) == 1) == 1); %The indices where the values start decreasing
startFallInd = decInds(decNonIncre(1));
end

function [endRiseInd, startFallInd] = findCrossingsManual(sensor, fs)

[B, A] = butter(8, (50/(fs/2)), 'low'); 
sensorFilt = filtfilt(B,A, sensor);

sensDiff = [0; diff(sensorFilt)]*20;
sensDiff2= [0; diff(sensDiff)]*20;

PresFig = figure;
plot(sensor, 'k'); hold on
plot(sensDiff, 'r'); hold on
plot(sensDiff2, 'g')
axis([0 3200 3.8 4.5])

[x1, ~] = getpts(PresFig);
[x2, ~] = getpts(PresFig);

endRiseInd = round(x1);
startFallInd = round(x2);

close all;
end

function lims = identifyLimits(niAn)

%Full Inidividual Trials: Pressure Sensor
lims.pressure    = [0 4 0 5];

%Aligned Pressure Data
lims.pressureAl = [0 3.5 0 5];

%Full Individual Trials: Force Sensors
lims.force       = [0 4 1 5];

%Full Individual Trials: f0 Audio 
lims.audio     = [0 4 -100 100];

%Section Mean Trials: f0 Audio 
lims.audioMean = [0 4 -100 100];

%Full trial f0 analysis

% %Sectioned f0 Analysis
% [~, Imax] = max(niAn.meanTrialf0_St(:,1,2)); %Mean Microphone f0, Perturbed Trials
% upBound_St = round(niAn.meanTrialf0_St(Imax,1,2) + niAn.meanTrialf0_St(Imax,2,2) + 10);
% [~, Imin] = min(niAn.meanTrialf0_St(:,1,2)); %Mean Microphone f0, Perturbed Trials
% lwBound_St = round(niAn.meanTrialf0_St(Imin,1,2) - niAn.meanTrialf0_St(Imin,2,2) - 10);
% 
% [~, Imax] = max(niAn.meanTrialf0_Sp(:,1,2)); %Mean Microphone f0, Perturbed Trials
% upBound_Sp = round(niAn.meanTrialf0_Sp(Imax,1,2) + niAn.meanTrialf0_Sp(Imax,2,2) + 10);
% [~, Imin] = min(niAn.meanTrialf0_Sp(:,1,2)); %Mean Microphone f0, Perturbed Trials
% lwBound_Sp = round(niAn.meanTrialf0_Sp(Imin,1,2) - niAn.meanTrialf0_Sp(Imin,2,2) - 10);
% 
% if upBound_St > upBound_Sp
%     lims.upBoundSec = upBound_St;
% else
%     lims.upBoundSec = upBound_Sp;
% end
% 
% if lwBound_St < lwBound_Sp
%     lims.lwBoundSec = lwBound_St;
% else
%     lims.lwBoundSec = lwBound_Sp;
% end
end

function res = packResults(niAn, lims)

res.subject = niAn.subject;
res.run     = niAn.run;
res.curSess = niAn.curSess;

res.numTrials     = niAn.numTrial;
res.numPertTrials = niAn.numPertTrials;
res.pertIdx       = niAn.pertIdx;
res.pertTrig      = niAn.pertTrig;

res.timeS      = niAn.time_DN;
res.sensorP    = niAn.sensorP_p; %Individual Processed perturbed trials. 
res.lagTimeP   = niAn.lagsPres;
res.lagTimePm  = niAn.meanLagTimeP;
res.riseTimeP  = niAn.riseTimeP;
res.riseTimePm = niAn.meanRiseTimeP;
res.rangeP     = niAn.rangePressures;
res.rangePm    = niAn.meanRangePressure;
res.limitsP    = lims.pressure;

res.timeSAl   = niAn.time_Al;
res.sensorPAl = niAn.sensorP_Al;
res.limitsPAl = lims.pressureAl;

res.timeA    = niAn.time_audio;
res.f0b      = niAn.f0b;

%Full Individual Trials: Mic/Head f0 Trace 
res.audioMf0TrialPert = niAn.audioMf0_p;
res.audioMf0TrialCont = niAn.audioMf0_c;
res.audioHf0TrialPert = niAn.audioHf0_p;
res.audioHf0TrialCont = niAn.audioHf0_c;
res.limitsA           = lims.audio;

%Sectioned Mean Trials: Mic/Head f0 Trace 
res.audioMf0MeanPert = niAn.audioMf0_meanp; % [MeanSig 90%CI]
res.audioMf0MeanCont = niAn.audioMf0_meanc;
res.audioHf0MeanPert = niAn.audioHf0_meanp;
res.audioHf0MeanCont = niAn.audioHf0_meanc;
res.limitsAmean      = lims.audioMean;
end