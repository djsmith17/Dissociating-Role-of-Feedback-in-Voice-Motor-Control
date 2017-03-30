function niAn = dfAnalysisNIDAQ(expParam, DAQin)

[r, c, n] = size(DAQin);
sRate = expParam.sRateQ;

niAn.numTrial = n;
niAn.numCh    = c;
niAn.dnSamp   = 10;
niAn.time     = 0:1/sRate:(r-1)/sRate;
niAn.pertSig  = squeeze(DAQin(:,1,:));
niAn.sensorFC = squeeze(DAQin(:,2,:));
niAn.sensorFN = squeeze(DAQin(:,3,:));
niAn.sensorP  = squeeze(DAQin(:,4,:));
niAn.audioM   = squeeze(DAQin(:,5,:));
niAn.audioH   = squeeze(DAQin(:,6,:));
niAn.sensorO  = squeeze(DAQin(:,7,:));

[B,A] = butter(4, 40/(sRate/2)); %Low-pass filter under 40
niAn.sensorFC = filter(B,A,abs(niAn.sensorFC));
niAn.sensorFN = filter(B,A,abs(niAn.sensorFN));

niAn.pertSig_DN  = dnSampleSignal(niAn.pertSig, niAn.dnSamp);
niAn.sensorFC_DN = dnSampleSignal(niAn.sensorFC, niAn.dnSamp);
niAn.sensorFN_DN = dnSampleSignal(niAn.sensorFN, niAn.dnSamp);
niAn.sensorP_DN  = dnSampleSignal(niAn.sensorP, niAn.dnSamp);

[pertTrig, pertThresh, idxPert] = findPertTrigs(niAn.time, niAn.pertSig);
[presTrig, presThresh, idxPres] = findPertTrigs(niAn.time, niAn.sensorP);
[fSCTrig, fSCThresh, idxFC]     = findPertTrigs(niAn.time, niAn.sensorFC);  
[fSNTrig, fSNThresh, idxFN]     = findPertTrigs(niAn.time, niAn.sensorFN); 

[niAn.lagsPres, niAn.lagMeansPres] = calcMeanLags(pertTrig, presTrig);
[niAn.lagsFC, niAn.lagMeansFC]     = calcMeanLags(pertTrig, fSCTrig);
[niAn.lagsFN, niAn.lagMeansFN]     = calcMeanLags(pertTrig, fSNTrig);

niAn.rangePressures = [];
for ii = 1:numTrial
    onsetPressure  = round(100*max(niAn.sensorP(:,ii)))/100;
    offsetPressure = round(100*niAn.sensorP(idxPert(ii,2), ii))/100;
    niAn.rangePressures = cat(1, niAn.rangePressures, [onsetPressure offsetPressure]);
end
niAn.meanRangePressure = mean(niAn.rangePressures, 1);

niAn.pSensorAl = alignSensorData(sRate, niAn.numTrial, niAn.pertidx, niAn.sensorP);
niAn.timeAl    = 0:1/sRate:(length(niAn.pSensorAl)-1)/sRate;

niAn.trigs      = pertTrig;
niAn.pertThresh = pertThresh;
niAn.presTrig   = presTrig;
niAn.presThresh = presThresh;
niAn.fSCTrig    = fSCTrig;
niAn.fSCThresh  = fSCThresh;
niAn.fSNTrig    = fSNTrig;
niAn.fSNThresh  = fSNThresh;
end

function [trigs, threshes, idx] = findPertTrigs(time, sensor)
sensor = round(sensor); %Should be step function 0V or 3V
[~, numTrial] = size(sensor);

trigs = [];
threshes = [];
idx   = [];
for i = 1:numTrial
    thresh = mean(sensor(2000:4000, i));
    
    I = find(sensor(:,i) > thresh);
    trigSt = round(1000*time(I(1)))/1000;
    trigSp = round(1000*time(I(end)))/1000;

    trigs    = cat(1, trigs, [trigSt trigSp]);
    threshes = cat(1, threshes, thresh);
    idx      = cat(1, idx, [I(1) I(end)]);
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