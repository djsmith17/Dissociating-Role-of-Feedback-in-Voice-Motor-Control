function niAn = dfAnalysisNIDAQ(expParam, DAQin)

[r, c, n] = size(DAQin);
sRate = expParam.sRateQ;

niAn.sRate    = sRate;
niAn.numTrial = n;
niAn.numCh    = c;
niAn.dnSamp   = 10;
niAn.sRateDN  = sRate/niAn.dnSamp;
niAn.time     = (0:1/sRate:(r-1)/sRate)';
niAn.pertSig  = squeeze(DAQin(:,1,:));
niAn.sensorFC = squeeze(DAQin(:,2,:));
niAn.sensorFN = squeeze(DAQin(:,3,:));
niAn.sensorP  = squeeze(DAQin(:,4,:));
niAn.audioM   = squeeze(DAQin(:,5,:));
niAn.audioH   = squeeze(DAQin(:,6,:));
niAn.sensorO  = squeeze(DAQin(:,7,:));

[B,A] = butter(4, 10/(sRate/2)); %Low-pass filter under 40
niAn.sensorFC = filter(B,A,abs(niAn.sensorFC));
niAn.sensorFN = filter(B,A,abs(niAn.sensorFN));

niAn.time_DN     = dnSampleSignal(niAn.time, niAn.dnSamp);
niAn.pertSig_DN  = dnSampleSignal(niAn.pertSig, niAn.dnSamp);
niAn.sensorP_DN  = dnSampleSignal(niAn.sensorP, niAn.dnSamp);
niAn.sensorFC_DN = dnSampleSignal(niAn.sensorFC, niAn.dnSamp);
niAn.sensorFN_DN = dnSampleSignal(niAn.sensorFN, niAn.dnSamp);

[niAn.pertTrig, niAn.idxPert] = findPertTrigs(niAn.time_DN, niAn.pertSig_DN, niAn.sRateDN);
[niAn.presTrig, niAn.idxPres] = findPertTrigs(niAn.time_DN, niAn.sensorP_DN, niAn.sRateDN);
[niAn.fSCTrig, niAn.idxFC]    = findPertTrigs(niAn.time_DN, niAn.sensorFC_DN, niAn.sRateDN);  
[niAn.fSNTrig, niAn.idxFN]    = findPertTrigs(niAn.time_DN, niAn.sensorFN_DN, niAn.sRateDN); 

[niAn.lagsPres, niAn.lagMeansPres] = calcMeanLags(niAn.pertTrig, niAn.presTrig);
[niAn.lagsFC, niAn.lagMeansFC]     = calcMeanLags(niAn.pertTrig, niAn.fSCTrig);
[niAn.lagsFN, niAn.lagMeansFN]     = calcMeanLags(niAn.pertTrig, niAn.fSNTrig);

niAn.rangePressures = [];
niAn.maxTimesP       = [];
for ii = 1:niAn.numTrial
    [maxP, maxInd] = max(niAn.sensorP_DN(:,ii));
    [minP  ] = niAn.sensorP_DN(niAn.idxPert(ii,2), ii);      
    
    onsetPressure  = round(100*maxP)/100;
    offsetPressure = round(100*minP)/100;
    maxTime = round(100*niAn.time_DN(maxInd))/100;
    niAn.maxTimesP      = cat(1, niAn.maxTimesP, maxTime);
    niAn.rangePressures = cat(1, niAn.rangePressures, [onsetPressure offsetPressure]);
end
niAn.meanRangePressure = mean(niAn.rangePressures, 1);

niAn.riseTimeP = niAn.maxTimesP - niAn.presTrig(:,1);
niAn.meanRiseTimeP = mean(niAn.riseTimeP);

niAn.sensorP_Al = alignSensorData(niAn.sRateDN , niAn.numTrial, niAn.idxPert, niAn.sensorP_DN);
niAn.time_Al    = 0:1/niAn.sRateDN :(length(niAn.sensorP_Al)-1)/niAn.sRateDN ;
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

function sensorDN = dnSampleSignal(sensor, dnSamp)
[numSamp, numTrial] = size(sensor);
numSampDN = numSamp/dnSamp;

sensorDN = zeros(numSampDN, numTrial);
for i = 1:numSampDN
    sensorDN(i,:) = mean(sensor((1:dnSamp) + dnSamp*(i-1),:));
end
end

function sensorDN = dnSampleSmoothSignal(sensor, fs)
[numSamp, ~] = size(sensor);

win = fs*0.075;
pOV = 0.8;
tStep = win*(1-pOV);
stIdx = 1:tStep:numSamp-win;

sensorDN = [];
for iSt = stIdx
    sensorDN = cat(1, sensorDN, mean(sensor(iSt:iSt+win, :)));
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

function sensorAl = alignSensorData(sRate, numTrial, idx, sensor)

sensorAl = [];
for ii = 1:numTrial
    St = idx(ii,1) - sRate*1;
    Sp = idx(ii,1) + sRate*2.5;
    
    sensorAl = cat(2, sensorAl, sensor(St:Sp,ii));
end
end