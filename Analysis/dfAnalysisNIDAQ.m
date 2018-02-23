function [niAn, niRes] = dfAnalysisNIDAQ(dirs, expParam, DAQin, bTf0b, AudFlag, PresFlag, iRF)
%A quick reference
%
%Pert: Perturbation signal
%P:    Pressure sensor signal
%FC:   Force Sensor Collar
%FN:   Force Sensor Neck
%
%Trig: Trigger values where onset and offset occur
%DN:   Down Sampled (and smoothed)

sv2File = 0;
[r, c, n] = size(DAQin);
sRate = expParam.sRateQ;

%Identify some starting variables
niAn.AnaType  = 'NIDAQ';
niAn.expType  = expParam.expType;
niAn.subject  = expParam.subject;
niAn.run      = expParam.run;
niAn.curSess  = expParam.curSess;
niAn.f0AnaFile = [niAn.subject niAn.run 'f0AnalysisN.mat'];
niAn.gender   = expParam.gender;
niAn.AudFB    = expParam.AudFB;
niAn.bTf0b    = bTf0b;
niAn.trialType = expParam.trialType;

fprintf('Starting NIDAQ Analysis for %s, %s with f0 of %0.2f Hz\n', niAn.subject, niAn.run, niAn.bTf0b)

niAn.dnSamp   = 10;
niAn.sRate    = sRate;
niAn.numSamp  = r;
niAn.numTrial = n;
niAn.numCh    = c;
niAn.expTrigs = expParam.trigs(:,:,1); %Time

%Find all the perturbed trials
[niAn.ContTrials, niAn.contIdx] = find(niAn.trialType == 0);
[niAn.PertTrials, niAn.pertIdx] = find(niAn.trialType == 1);
niAn.numContTrials = sum(niAn.ContTrials);
niAn.numPertTrials = sum(niAn.PertTrials);

%Unpack the NIDAQ raw data set
niAn.time     = (0:1/niAn.sRate:(niAn.numSamp-1)/niAn.sRate)';
niAn.pertSig  = squeeze(DAQin(:,1,:));
niAn.sensorFC = squeeze(DAQin(:,2,:));
niAn.sensorFN = squeeze(DAQin(:,3,:));
niAn.sensorP  = squeeze(DAQin(:,4,:));
niAn.audioM   = squeeze(DAQin(:,5,:));
niAn.audioH   = squeeze(DAQin(:,6,:));
niAn.sensorO  = squeeze(DAQin(:,7,:));

%ZeroMean the Offset
niAn.sensorPz = correctBaseline(niAn.sensorP, niAn.sRate);

%Preprocessing some of the Force sensors
niAn.sensorFCz = sensorPreProcessing(niAn.sensorFC, niAn.sRate);
niAn.sensorFNz = sensorPreProcessing(niAn.sensorFN, niAn.sRate);

niAn.sRateDN     = niAn.sRate/niAn.dnSamp;
niAn.time_DN     = dnSampleSignal(niAn.time, niAn.dnSamp);    % DownSampled Time
niAn.pertSig_DN  = dnSampleSignal(niAn.pertSig, niAn.dnSamp); % DownSampled Perturbatron Signal
niAn.sensorP_DN  = dnSampleSignal(niAn.sensorPz, niAn.dnSamp);
niAn.sensorFC_DN = dnSampleSignal(niAn.sensorFCz, niAn.dnSamp);
niAn.sensorFN_DN = dnSampleSignal(niAn.sensorFNz, niAn.dnSamp);

%Parse out the perturbed trials
niAn.pertSig_p  = parseTrialTypes(niAn.pertSig_DN, niAn.pertIdx);  % Only Perturbed Trials
niAn.sensorP_p  = parseTrialTypes(niAn.sensorP_DN, niAn.pertIdx); % Only Perturbed Trials
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

niAn.OnOfValP  = [];
niAn.OnOfValPm = [];
niAn.riseTimeP = [];
niAn.riseTimePm = [];
niAn.sensorP_Al = [];
niAn.time_Al    = [];
if PresFlag == 1
    %Sensor Dynamics of the Pressure Sensor
    [niAn.OnOfValP,  niAn.OnOfValPm, ...
     niAn.riseTimeP, niAn.riseTimePm] = ...
    analyzeSensorDynamics(niAn.time_DN, niAn.sensorP_p, niAn.sRateDN, niAn.presTrig);

    %Aligning pressure signal for perturbed trials
    niAn.sensorP_Al = alignSensorData(niAn.sensorP_p, niAn.sRateDN, niAn.idxPert);
    niAn.time_Al    = (0:1/niAn.sRateDN :(length(niAn.sensorP_Al)-1)/niAn.sRateDN)';
end

%The Audio Analysis
niAn = dfAnalysisAudio(dirs, niAn, AudFlag, iRF);
    
lims  = identifyLimits(niAn);
niRes = packResults(niAn, lims);
end

function sensorZeroed = correctBaseline(sensor, fs)
%The first second of the first trial will have the baseline pressure reading
%from the sensor. We will use this to fix the offset in all the trials. 

firstS  = 1:(1*fs);          % Grab the first second
firstT  = sensor(firstS, 1); % Grab the very first trial
meanRec = mean(firstT);      % I really mean it

sensorZeroed = sensor - meanRec;
end

function sensorPP     = sensorPreProcessing(sensor, sRate)
%This was mostly to mess around with the force sensor, but for right now we
%will hide that all in here. Likely will never need this again. 

sensorRed = 4*(sensor-2);

[B,A] = butter(4, 10/(sRate/2)); %Low-pass filter under 10
sensorPP = filter(B,A,abs(sensorRed));
end

function sensorDN = dnSampleSignal(sensor, dnSamp)
[numSamp, numTrial] = size(sensor);
numSampDN = numSamp/dnSamp;

sensorDN = zeros(numSampDN, numTrial);
for i = 1:numSampDN
    sensorDN(i,:) = mean(sensor((1:dnSamp) + dnSamp*(i-1),:));
end
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

sensDiff  = smooth([0; diff(sensorFilt)]*50, 10);

diffPeakMax = 0.05;
peakLeadLag = fs*0.5; % How many seconds past he peak to look for the level off. 

if man == 0
    
    %yes its a bit lazy. I am sorry. 
    [pksU, locU] = findpeaks(sensDiff);
    [pksD, locD] = findpeaks(-1*sensDiff);    

    if isempty(pksU)
        disp('No rise found')
        endRiseInd = [];
    else
        [~, maxInd] = max(pksU);
        maxIndFull = locU(maxInd);
        searchR = maxIndFull:maxIndFull+peakLeadLag;
        diffDownRamp = find(sensDiff(searchR) > diffPeakMax);
        
        endRiseInd = searchR(diffDownRamp(end));
    end
    
    if isempty(pksD)
        disp('No fall found')
        startFallInd = [];
    else
        [~, minInd] = max(pksD);
        minIndFull = locD(minInd);
        searchR = minIndFull-peakLeadLag:minIndFull;
        diffUpRamp = find(sensDiff(searchR) > diffPeakMax);
        
        startFallInd = searchR(diffUpRamp(1));
    end
    
%     figure
%     plot(sensor,'k')
%     hold on
%     plot([endRiseInd endRiseInd], [-100 100], 'b')
%     hold on
%     plot([startFallInd startFallInd], [-100 100], 'r')
%     axis([0 3200 -0.05 4])
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

function lims = identifyLimits(An)

%Full Inidividual Trials: Pressure Sensor
lims.pressure   = [0 4 0 5];

%Aligned Pressure Data
lims.pressureAl = [0 3.5 -0.5 5];

%Full Individual Trials: Force Sensors
lims.force      = [0 4 1 5];

%%%%%%%%%%%lims.audio%%%%%%%%%%%
%Full Individual Trials: f0 Audio
if ~isempty(An.audioMf0_pPP)
    pertTrialsM = An.audioMf0_pPP;
    pertTrialsH = An.audioHf0_pPP;
    sec = 100:700;

    uLMa = max(pertTrialsM(sec,:));
    uLMa(find(uLMa > 150)) = 0;
    lLMa = min(pertTrialsM(sec,:));
    lLMa(find(lLMa < -150)) = 0;

    uLM = round(max(uLMa)) + 20;
    lLM = round(min(lLMa)) - 20;
       
    uLHa = max(pertTrialsH(sec,:));
    uLHa(find(uLHa > 150)) = 0;
    lLHa = min(pertTrialsH(sec,:));
    lLHa(find(lLHa < -150)) = 0;
    
    uLH = round(max(uLHa)) + 20;
    lLH = round(min(lLHa)) - 20;
    
    if uLH > uLM
        uLMH = uLH;
    else
        uLMH = uLM;
    end

    if lLH < lLM
        lLMH = lLH;
    else
        lLMH = lLM;
    end
    
    lims.audioM         = [0 4 lLM uLM];
    lims.audioAudRespMH = [0 4 lLMH uLMH];
else
    lims.audioM         = [0 4 -20 20];
    lims.audioAudRespMH = [0 4 -100 100];
end

%Section Mean Pertrubed Trials: f0 Audio 
if ~isempty(An.audioMf0_meanp)
    [~, Imax] = max(An.audioMf0_meanp(:,1)); %Max Pert Onset
    upBoundOn = round(An.audioMf0_meanp(Imax,1) + An.audioMf0_meanp(Imax,2) + 10);
    [~, Imin] = min(An.audioMf0_meanp(:,1)); %Min Pert Onset
    lwBoundOn = round(An.audioMf0_meanp(Imin,1) - An.audioMf0_meanp(Imin,2) - 10);

    [~, Imax] = max(An.audioMf0_meanp(:,3)); %Max Pert Offset
    upBoundOf = round(An.audioMf0_meanp(Imax,3) + An.audioMf0_meanp(Imax,4) + 10);
    [~, Imin] = min(An.audioMf0_meanp(:,3)); %Min Pert Offset
    lwBoundOf = round(An.audioMf0_meanp(Imin,3) - An.audioMf0_meanp(Imin,4) - 10);

    if upBoundOn > upBoundOf
        upBoundM = upBoundOn;
    else
        upBoundM = upBoundOf;
    end

    if lwBoundOn < lwBoundOf
        lwBoundM = lwBoundOn;
    else
        lwBoundM = lwBoundOf;
    end   

    lims.audioMean = [-0.5 1.0 lwBoundM upBoundM];
else
    lims.audioMean = [-0.5 1.0 -50 50];
end

%%%%%%%%%%%lims.audioMH%%%%%%%%%%%%
%Section Mean Pertrubed Trials: f0 Audio 
if ~isempty(An.audioHf0_meanp)
    [~, Imax] = max(An.audioHf0_meanp(:,1)); %Max Pert Onset
    upBoundOn = round(An.audioHf0_meanp(Imax,1) + An.audioHf0_meanp(Imax,2) + 10);
    [~, Imin] = min(An.audioHf0_meanp(:,1)); %Min Pert Onset
    lwBoundOn = round(An.audioHf0_meanp(Imin,1) - An.audioHf0_meanp(Imin,2) - 10);

    [~, Imax] = max(An.audioHf0_meanp(:,3)); %Max Pert Offset
    upBoundOf = round(An.audioHf0_meanp(Imax,3) + An.audioHf0_meanp(Imax,4) + 10);
    [~, Imin] = min(An.audioHf0_meanp(:,3)); %Min Pert Offset
    lwBoundOf = round(An.audioHf0_meanp(Imin,3) - An.audioHf0_meanp(Imin,4) - 10);

    if upBoundOn > upBoundOf
        upBoundH = upBoundOn;
    else
        upBoundH = upBoundOf;
    end

    if lwBoundOn < lwBoundOf
        lwBoundH = lwBoundOn;
    else
        lwBoundH = lwBoundOf;
    end
    
    if upBoundH > upBoundM
        upBoundMH = upBoundH;
    else
        upBoundMH = upBoundM;
    end
    
    if lwBoundH < lwBoundM
        lwBoundMH = lwBoundH;
    else
        lwBoundMH = lwBoundM;
    end

    lims.audioMH = [-0.5 1.0 lwBoundMH upBoundMH];
else
    lims.audioMH = [-0.5 1.0 -100 50];
end

end

function res = packResults(niAn, lims)

res.expType = niAn.expType;
res.subject = niAn.subject;
res.run     = niAn.run;
res.curSess = niAn.curSess;
res.AudFB   = niAn.AudFB;

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

res.timef0    = niAn.timef0;
res.f0b       = niAn.f0b;

res.numContTrialsPP = niAn.numContTrialsPP;
res.numPertTrialsPP = niAn.numPertTrialsPP;
res.pertTrigPP      = niAn.pertTrigPP;

%Full Individual Trials: Mic/Head f0 Trace 
res.audioMf0TrialPert = niAn.audioMf0_pPP;
res.audioMf0TrialCont = niAn.audioMf0_cPP;
res.audioHf0TrialPert = niAn.audioHf0_pPP;
res.audioHf0TrialCont = niAn.audioHf0_cPP;
res.limitsA           = lims.audioM;
res.limitsAudRes      = lims.audioAudRespMH;

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
res.limitsAMH        = lims.audioMH;

%Inflation Response
res.respVar      = niAn.respVar;
res.respVarM     = niAn.respVarM;
res.respVarSD    = niAn.respVarSD;
res.InflaStimVar = niAn.InflaStimVar;
end