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

niAn.pertSig_DN  = dnSampleSignal(niAn.pertSig, niAn.dnSamp, sRate);
niAn.sensorFC_DN = dnSampleSignal(niAn.sensorFC, niAn.dnSamp, sRate);
niAn.sensorFN_DN = dnSampleSignal(niAn.sensorFN, niAn.dnSamp, sRate);
niAn.sensorP_DN  = dnSampleSignal(niAn.sensorP, niAn.dnSamp, sRate);

[niAn.pertTrig, niAn.pertThresh, niAn.idxPert] = findPertTrigs(niAn.time, niAn.pertSig_DN, sRate);
[niAn.presTrig, niAn.presThresh, niAn.idxPres] = findPertTrigs(niAn.time, niAn.sensorP_DN, sRate);
[niAn.fSCTrig, niAn.fSCThresh, niAn.idxFC]     = findPertTrigs(niAn.time, niAn.sensorFC_DN, sRate);  
[niAn.fSNTrig, niAn.fSNThresh, niAn.idxFN]     = findPertTrigs(niAn.time, niAn.sensorFN_DN, sRate); 

[niAn.lagsPres, niAn.lagMeansPres] = calcMeanLags(niAn.pertTrig, niAn.presTrig);
[niAn.lagsFC, niAn.lagMeansFC]     = calcMeanLags(niAn.pertTrig, niAn.fSCTrig);
[niAn.lagsFN, niAn.lagMeansFN]     = calcMeanLags(niAn.pertTrig, niAn.fSNTrig);

niAn.rangePressures = [];
for ii = 1:numTrial
    onsetPressure  = round(100*max(niAn.sensorP(:,ii)))/100;
    offsetPressure = round(100*niAn.sensorP(niAn.idxPert(ii,2), ii))/100;
    niAn.rangePressures = cat(1, niAn.rangePressures, [onsetPressure offsetPressure]);
end
niAn.meanRangePressure = mean(niAn.rangePressures, 1);

niAn.pSensorAl = alignSensorData(sRate, niAn.numTrial, niAn.pertidx, niAn.sensorP);
niAn.timeAl    = 0:1/sRate:(length(niAn.pSensorAl)-1)/sRate;
end

function [trigs, threshes, idx] = findPertTrigs(time, sensor, fs)
[numSamp, numTrial] = size(sensor);
st = 1*fs;
sp = 3*fs;

trigs = [];
threshes = [];
idx   = [];
for i = 1:numTrial
    for j = 2:numSamp
        m(j) = sensor(j,i) - sensor((j-1),i);
    end
        
    thresh = mean(sensor(1:st, i));
    
    I = find(sensor(st:sp,i) > thresh);
    trigSt = round(1000*time(I(1)))/1000;
    trigSp = round(1000*time(I(end)))/1000;

    trigs    = cat(1, trigs, [trigSt trigSp]);
    threshes = cat(1, threshes, thresh);
    idx      = cat(1, idx, [I(1) I(end)]);
end
end

function sensorDN = dnSampleSignal(sensor, dnSamp, fs)
[numSamp, numTrial] = size(sensor);

win = fs*0.075;
numSampDN = numSamp/dnSamp;

sensorDN = zeros(numSampDN, numTrial);
for i = 1:numSampDN
    sensorDN(i,:) = mean(sensor((1:win) + dnSamp*(i-1),:));
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