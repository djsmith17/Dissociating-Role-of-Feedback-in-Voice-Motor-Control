function niAn = dfAnalysisNIDAQ(expParam, DAQin)

[r, c, n] = size(DAQin);
sRate = expParam.sRateQ;

niAn.numTrial = n;
niAn.numCh    = c;
niAn.time     = 0:1/sRate:(r-1)/sRate;
niAn.pert     = squeeze(DAQin(:,1,:));
niAn.fSensorC = squeeze(DAQin(:,2,:));
niAn.fSensorN = squeeze(DAQin(:,3,:));
niAn.pSensor  = squeeze(DAQin(:,4,:));

[B,A] = butter(4, 40/(sRate/2)); %Low-pass filter under 40
niAn.fSensorC  = filter(B,A,abs(niAn.fSensorC));
niAn.fSensorN  = filter(B,A,abs(niAn.fSensorN));

[pertTrig, pertThresh, pertidx] = findPertTrigs(niAn.time, niAn.pert);
[presTrig, presThresh, presidx] = findPertTrigs(niAn.time, niAn.pSensor);
[fSCTrig, fSCThresh, fSCidx]    = findPertTrigs(niAn.time, niAn.fSensorC);  
[fSNTrig, fSNThresh, fSNidx]    = findPertTrigs(niAn.time, niAn.fSensorN); 

[PresLags, PresLagVals] = calcMeanLags(pertTrig, presTrig);
[fSCLags, fSCLagVals]   = calcMeanLags(pertTrig, fSCTrig);
[fSNLags, fSNLagVals]   = calcMeanLags(pertTrig, fSNTrig);

rangePressures = [];
for ii = 1:numTrial
    onsetPressure  = round(100*max(niAn.pSensor(:,ii)))/100;
    offsetPressure = round(100*niAn.pSensor(pertidx(ii,2), ii))/100;
    rangePressures = cat(1, rangePressures, [onsetPressure offsetPressure]);
end

niAn.pSensorAl = alignSensorData(sRate, niAn.numTrial, niAn.pertidx, niAn.pSensor);
niAn.timeAl    = 0:1/sRate:(length(niAn.pSensorAl)-1)/sRate;

niAn.trigs      = pertTrig;
niAn.pertThresh = pertThresh;
niAn.presTrig   = presTrig;
niAn.presThresh = presThresh;
niAn.fSCTrig    = fSCTrig;
niAn.fSCThresh  = fSCThresh;
niAn.fSNTrig    = fSNTrig;
niAn.fSNThresh  = fSNThresh;

niAn.PresLags    = PresLags;
niAn.PresLagVals = PresLagVals;
niAn.fSCLags     = fSCLags;
niAn.fSCLagVals  = fSCLagVals;
niAn.fSNLags     = fSNLags;
niAn.fSNLagVals  = fSNLagVals;

niAn.rangePressures = rangePressures;
niAn.meanRangePressure = mean(rangePressures, 1);
end

function [trigs, threshes, idx] = findPertTrigs(time, pertCh)
pertCh = round(pertCh); %Should be step function 0V or 3V
[~, c] = size(pertCh);

trigs = [];
threshes = [];
idx   = [];
for i = 1:c
    thresh = mean(pertCh(2000:4000, i));
    
    I = find(pertCh(:,i) > thresh);
    trigSt = round(1000*time(I(1)))/1000;
    trigSp = round(1000*time(I(end)))/1000;

    trigs    = cat(1, trigs, [trigSt trigSp]);
    threshes = cat(1, threshes, thresh);
    idx      = cat(1, idx, [I(1) I(end)]);
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