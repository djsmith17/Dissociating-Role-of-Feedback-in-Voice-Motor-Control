function [niAn, res] = dfAnalysisNIDAQ(dirs, expParam, DAQin, audioFlag)
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
[niAn.ContTrials, niAn.contIdx] = find(niAn.trialType == 0);
[niAn.PertTrials, niAn.pertIdx] = find(niAn.trialType == 1);
niAn.numContTrials = sum(niAn.ContTrials);
niAn.numPertTrials = sum(niAn.PertTrials);

%Identify a few analysis varaibles
fV.win    = 0.05;  %seconds
fV.fsA    = 1/fV.win;
fV.winP   = fV.win*niAn.sRate;
fV.pOV    = 0.60;  %60% overlap
fV.tStepP = fV.winP*(1-fV.pOV);
fV.winSts = 1:fV.tStepP:(niAn.numSamp-fV.winP);
fV.numWin = length(fV.winSts);
fV.freqCutOff = 300;
niAn.fV = fV;

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
niAn.pertSig_DN  = dnSampleSignal(niAn.pertSig, niAn.dnSamp); % DownSampled Perturbatron Signal
niAn.sensorP_DN  = dnSampleSignal(niAn.sensorP, niAn.dnSamp);
niAn.sensorFC_DN = dnSampleSignal(niAn.sensorFC, niAn.dnSamp);
niAn.sensorFN_DN = dnSampleSignal(niAn.sensorFN, niAn.dnSamp);

%ZeroMean the Offset
niAn.sensorP_DNz = correctPOffset(niAn.sensorP_DN, niAn.sRateDN);

%Parse out the pertrubed trials
niAn.pertSig_p  = parseTrialTypes(niAn.pertSig_DN, niAn.pertIdx);  % Only Perturbed Trials
niAn.sensorP_p  = parseTrialTypes(niAn.sensorP_DNz, niAn.pertIdx); % Only Perturbed Trials
niAn.sensorFC_p = parseTrialTypes(niAn.sensorFC_DN, niAn.pertIdx); % Only Perturbed Trials
niAn.sensorFN_p = parseTrialTypes(niAn.sensorFN_DN, niAn.pertIdx); % Only Perturbed Trials

%Make a dummy set of contTrig
niAn.contTrig = repmat([1 2.5], niAn.numContTrials, 1);

%Find Rising and Falling Edges of sensor signals
[niAn.pertTrig, niAn.idxPert] = findPertTrigs(niAn.time_DN, niAn.pertSig_p, niAn.sRateDN);
[niAn.presTrig, niAn.idxPres] = findPertTrigs(niAn.time_DN, niAn.sensorP_p, niAn.sRateDN);
[niAn.fSCTrig, niAn.idxFC]    = findPertTrigs(niAn.time_DN, niAn.sensorFC_p, niAn.sRateDN);  
[niAn.fSNTrig, niAn.idxFN]    = findPertTrigs(niAn.time_DN, niAn.sensorFN_p, niAn.sRateDN); 

[niAn.lagsPres, niAn.meanLagTimeP] = calcMeanLags(niAn.pertTrig, niAn.presTrig);
[niAn.lagsFC, niAn.meanLagTimeFC]  = calcMeanLags(niAn.pertTrig, niAn.fSCTrig);
[niAn.lagsFN, niAn.meanLagTimeFN]  = calcMeanLags(niAn.pertTrig, niAn.fSNTrig);

%Sensor Dynamics of the Pressure Sensor
[niAn.OnOfValP,  niAn.OnOfValPm, ...
 niAn.riseTimeP, niAn.riseTimePm] = ...
analyzeSensorDynamics(niAn.time_DN, niAn.sensorP_p, niAn.sRateDN, niAn.presTrig);

%Aligning pressure signal for perturbed trials
niAn.sensorP_Al = alignSensorData(niAn.sensorP_p, niAn.sRateDN, niAn.idxPert);
niAn.time_Al    = (0:1/niAn.sRateDN :(length(niAn.sensorP_Al)-1)/niAn.sRateDN)';

if audioFlag == 1
    % Audio Processing
    
    [niAn.time_audio, niAn.audioMf0, niAn.fsA] = signalFrequencyAnalysis(dirs, niAn.time, niAn.audioM, niAn.sRate, fV, 1);
    [niAn.time_audio, niAn.audioHf0, niAn.fsA] = signalFrequencyAnalysis(dirs, niAn.time, niAn.audioH, niAn.sRate, fV, 1);
    prePert         = (0.5 < niAn.time_audio & 1.0 > niAn.time_audio);
    niAn.trialf0b   = mean(niAn.audioMf0(prePert,:),1);
    niAn.f0b        = mean(niAn.trialf0b);

    niAn.audioMf0_norm = normalizeDAQf0(niAn.audioMf0, niAn.trialf0b);
    niAn.audioHf0_norm = normalizeDAQf0(niAn.audioHf0, niAn.trialf0b);
    %Find the Perturbed Trials
    niAn.audioMf0_p = parseTrialTypes(niAn.audioMf0_norm, niAn.pertIdx);
    niAn.audioHf0_p = parseTrialTypes(niAn.audioHf0_norm, niAn.pertIdx);
    niAn.audioMf0_c = parseTrialTypes(niAn.audioMf0_norm, niAn.contIdx);
    niAn.audioHf0_c = parseTrialTypes(niAn.audioHf0_norm, niAn.contIdx);
    
    [niAn.secTime, niAn.audioMf0_Secp] = sectionAudioData(niAn.time_audio, niAn.audioMf0_p, niAn.fsA, niAn.pertTrig);
    [niAn.secTime, niAn.audioHf0_Secp] = sectionAudioData(niAn.time_audio, niAn.audioHf0_p, niAn.fsA, niAn.pertTrig);
    [niAn.secTime, niAn.audioMf0_Secc] = sectionAudioData(niAn.time_audio, niAn.audioMf0_c, niAn.fsA, niAn.contTrig);
    [niAn.secTime, niAn.audioHf0_Secc] = sectionAudioData(niAn.time_audio, niAn.audioHf0_c, niAn.fsA, niAn.contTrig);
    
    %Mean around the onset and offset
    niAn.audioMf0_meanp = meanAudioData(niAn.audioMf0_Secp);
    niAn.audioHf0_meanp = meanAudioData(niAn.audioHf0_Secp);
    niAn.audioMf0_meanc = meanAudioData(niAn.audioMf0_Secc);
    niAn.audioHf0_meanc = meanAudioData(niAn.audioHf0_Secc); 
    
    %The Inflation Response
    [niAn.respVar, niAn.respVarMean] = InflationResponse(niAn.secTime, niAn.audioMf0_Secp);
    
else
    niAn.time_audio     = []; niAn.fsA            = [];
    niAn.audioMf0       = []; niAn.audioHf0       = [];
    niAn.trialf0b       = []; niAn.f0b            = [];
    niAn.audioMf0_norm  = []; niAn.audioHf0_norm  = [];
    niAn.audioMf0_p     = []; niAn.audioHf0_p     = [];
    niAn.audioMf0_c     = []; niAn.audioHf0_c     = [];
    niAn.secTime        = [];
    niAn.audioMf0_Secp  = []; niAn.audioHf0_Secp  = [];
    niAn.audioMf0_Secc  = []; niAn.audioHf0_Secc  = [];
    niAn.audioMf0_meanp = []; niAn.audioHf0_meanp = [];
    niAn.audioMf0_meanc = []; niAn.audioHf0_meanc = [];
    niAn.respVar        = []; niAn.respVarMean    = [];
end
    
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

function sensorZeroed = correctPOffset(sensor, fs)
%The first second of the trial will have the baseline pressure reading
%from the sensor. We will use this to fix the offset in all the trials. 

firstS  = 1:(1*fs);
firstT  = sensor(firstS, 1);
meanRec = mean(firstT);

sensorZeroed = sensor - meanRec;
end

function signalParse = parseTrialTypes(signal, idx)
%Expects trials to be in columns 

signalParse = signal(:, idx); %This is a little lazy I know. Get over it. 
end

function [trigs, idx] = findPertTrigs(time, sensor, fs)
%findPertTrigs(time, sensor, fs) finds rising and falling edges in sensor
%data. It is expected that these signals will be mostly step functions
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

function [OnOfVals, OnOfValsMean, riseTime, riseTimeMean] = analyzeSensorDynamics(time, sensor, fs, sensTrig)
%Analyzing the dynamics of the sensor during onset and offset of the
%stimulus signal.
[~, numTrial] = size(sensor);

OnOfInd   = [];
OnOfTime  = [];
OnOfVals  = [];
for ii = 1:numTrial
    [endRiseInd, startFallInd] = findCrossings(sensor(:,ii), fs, 0);

    onsetTime      = round(100*time(endRiseInd))/100;
    offsetTime     = round(100*time(startFallInd))/100;
    
    onsetSensorVal  = round(100*sensor(endRiseInd, ii))/100;
    offsetSensorVal = round(100*sensor(startFallInd, ii))/100;
  
    OnOfInd   = cat(1, OnOfInd, [endRiseInd, startFallInd]);
    OnOfTime  = cat(1, OnOfTime, [onsetTime offsetTime]);
    OnOfVals  = cat(1, OnOfVals, [onsetSensorVal offsetSensorVal]);
end

OnOfValsMean = mean(OnOfVals, 1);
riseTime     = OnOfTime(:,1) - sensTrig(:,1);
riseTimeMean = mean(riseTime);
end

function [endRiseInd, startFallInd] = findCrossings(sensor, fs, man)

[B, A] = butter(8, (50/(fs/2)), 'low'); 
sensorFilt = filtfilt(B,A, sensor);

sensDiff  = [0; diff(sensorFilt)]*50;
sensDiff2 = [0; diff(sensDiff)]*50;

if man == 0
    %This is a bit lazy, I know. For right now it will work.
    rangeU = 1*fs:2*fs;
    rangeD = 2*fs:3*fs;
    
    %yes its a bit lazy. I am sorry. 
    [pksU, locU] = findpeaks(sensDiff2, 'MinPeakHeight', 6);
    [pksD, locD] = findpeaks(-1*sensDiff2, 'MinPeakHeight', 10);

    if isempty(pksU)
        disp('No rise found')
        endRiseInd = [];
    else
        endRiseInd = locU(end);
    end
    
    if isempty(pksD)
        disp('No fall found')
        startFallInd = [];
    else
        startFallInd = locD(1);
    end
    
else
    PresFig = figure;
    plot(sensor, 'k'); hold on
    plot(sensDiff, 'r'); hold on
    plot(sensDiff2, 'g')
    axis([0 3200 3.8 4.5])

    [x1, ~] = getpts(PresFig);
    [x2, ~] = getpts(PresFig);

    endRiseInd   = round(x1);
    startFallInd = round(x2);

    close PresFig;
end
end

function sensorAl = alignSensorData(sensor, fs, idx)
[~, numTrial] = size(sensor);

sensorAl = [];
for ii = 1:numTrial
    St = idx(ii,1) - fs*1;
    Sp = idx(ii,1) + fs*2.5;
    
    sensorAl = cat(2, sensorAl, sensor(St:Sp, ii));
end
end

function [timef0, audiof0, fsA] = signalFrequencyAnalysis(dirs, time, audio, fs, fV, flag)
[~, numTrial] = size(audio);

if flag == 1
    [timef0, audiof0, fsA] = dfCalcf0Praat(dirs, audio, fs);
else
    %Low-Pass filter for the given cut off frequency
    [B,A]    = butter(4,(fV.freqCutOff)/(fs/2));

    audiof0 = zeros(fV.numWin, numTrial);
    for j = 1:numTrial %Trial by Trial
        sensorHP = filtfilt(B,A,audio(:,j));
        
        timef0 = [];
        for i = 1:fV.numWin
            winIdx = fV.winSts(i):fV.winSts(i)+ fV.winP - 1;
            timef0 = cat(1, timef0, mean(time(winIdx)));
            f0 = dfCalcf0Chile(sensorHP(winIdx), fs);
            if isnan(f0)
%                 disp('hit')
                f0 = 100;
            end
            audiof0(i,j) = f0;
        end
    end
    fsA = fV.fsA;
end
end

function audio_norm = normalizeDAQf0(audio, f0b)
[~, numTrial] = size(audio);

audio_norm = [];
for ii = 1:numTrial
    audio_trial = 1200*log2(audio(:,ii)./f0b(ii));
    audio_norm  = cat(2, audio_norm, audio_trial);
end
end

function [secTime, secAudio] = sectionAudioData(time, audio, fs, trigs)
[~, numTrial] = size(audio);
preEve  = 0.5; posEve = 1.0;
per     = 1/fs;
preEveP = preEve*fs;
posEveP = posEve*fs-1;

trigsR   = round2matchfs(trigs);

secAudio   = [];
OnsetSecs  = [];
OffsetSecs = [];
for ii = 1:numTrial 
    OnsetI  = find(time == trigsR(ii,1));
    OffsetI = find(time == trigsR(ii,2));
    
    OnsetISt = OnsetI - preEveP;
    OnsetISp = OnsetI + posEveP;
    
    OffsetISt = OffsetI - preEveP;
    OffsetISp = OffsetI + posEveP;
        
    OnsetSec  = audio(OnsetISt:OnsetISp, ii);
    OffsetSec = audio(OffsetISt:OffsetISp, ii);
    
    OnsetSecs  = cat(2, OnsetSecs, OnsetSec);
    OffsetSecs = cat(2, OffsetSecs, OffsetSec);
end

secTime  = (-0.5+per):per:1.0;
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

function [respVar, respVarm] = InflationResponse(secTime, secAudio)
[~, numTrial, ~] = size(secAudio);
postOnset = find(0 <= secTime & .20 >= secTime); % Cheating

stimMag  = [];
respMag  = [];
respPer  = [];
stimMagT = [];
for i = 1:numTrial
    onset = secAudio(:,i,1);    
    [minOn, minIdx] = min(onset(postOnset)); 
    stimMag  = cat(1, stimMag, minOn);       % Cheating
    respMag  = cat(1, respMag, onset(end));  % Cheating
    respPer  = cat(1, respPer, 100*(1 - onset(end)/minOn));
    stimMagT = cat(1, stimMagT, secTime(postOnset(minIdx)));
end

respVar  = [stimMagT stimMag respMag respPer];
respVarm = mean(respVar, 1);
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
[~, Imax] = max(niAn.audioMf0_meanp(:,1)); %Max Pert Onset
upBoundOn = round(niAn.audioMf0_meanp(Imax,1) + niAn.audioMf0_meanp(Imax,2) + 10);
[~, Imin] = min(niAn.audioMf0_meanp(:,1)); %Min Pert Onset
lwBoundOn = round(niAn.audioMf0_meanp(Imin,1) - niAn.audioMf0_meanp(Imin,2) - 10);

[~, Imax] = max(niAn.audioMf0_meanp(:,3)); %Max Pert Offset
upBoundOf = round(niAn.audioMf0_meanp(Imax,3) + niAn.audioMf0_meanp(Imax,4) + 10);
[~, Imin] = min(niAn.audioMf0_meanp(:,3)); %Min Pert Offset
lwBoundOf = round(niAn.audioMf0_meanp(Imin,3) - niAn.audioMf0_meanp(Imin,4) - 10);

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

function res = packResults(niAn, lims)

res.subject = niAn.subject;
res.run     = niAn.run;
res.curSess = niAn.curSess;

res.numTrials     = niAn.numTrial;
res.numContTrials = niAn.numContTrials;
res.numPertTrials = niAn.numPertTrials;
res.contIdx       = niAn.contIdx;
res.pertIdx       = niAn.pertIdx;
res.pertTrig      = niAn.pertTrig;

res.timeS      = niAn.time_DN;
res.sensorP    = niAn.sensorP_p; %Individual Processed perturbed trials. 
res.lagTimeP   = niAn.lagsPres;
res.lagTimePm  = niAn.meanLagTimeP;
res.riseTimeP  = niAn.riseTimeP;
res.riseTimePm = niAn.riseTimePm;
res.OnOfValP   = niAn.OnOfValP;
res.OnOfValPm  = niAn.OnOfValPm;
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

%Sections Trials: Mic/Head f0
res.secTime          = niAn.secTime;
res.audioMf0SecPert  = niAn.audioMf0_Secp;
res.audioMf0SecCont  = niAn.audioMf0_Secc;
res.audioHf0SecPert  = niAn.audioHf0_Secp;
res.audioHf0SecCont  = niAn.audioHf0_Secc;

%Mean Sectioned Trials: Mic/Head f0 Trace 
res.audioMf0MeanPert = niAn.audioMf0_meanp; % [MeanSigOn 90%CI MeanSigOff 90%CI]
res.audioMf0MeanCont = niAn.audioMf0_meanc;
res.audioHf0MeanPert = niAn.audioHf0_meanp;
res.audioHf0MeanCont = niAn.audioHf0_meanc;
res.limitsAmean      = lims.audioMean;

%Inflation Response
res.respVar  = niAn.respVar;
res.respVarM = niAn.respVarMean;
end