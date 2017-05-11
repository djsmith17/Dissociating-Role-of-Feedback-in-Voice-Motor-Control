function niAn = dfAnalysisNIDAQ(expParam, DAQin)
%A quick reference
%
%Pert: Perturbation signal
%P:    Pressure sensor signal
%FC:   Force Sensor Collar
%FN:   Force Sensor Neck
%
%Trig: Trigger values where onset and offset occur
%DN:   Down Sampled (and smoothed)

[r, c, n] = size(DAQin);
sRate = expParam.sRateQ;

niAn.subject  = expParam.subject;
niAn.run      = expParam.run;
niAn.curExp   = expParam.curExp;
niAn.sRate    = sRate;
niAn.numSamp  = r;
niAn.numTrial = n;
niAn.numCh    = c;
niAn.dnSamp   = 10;

niAn.win = 0.05;  %seconds
niAn.pOV = 0.60;  %60% overlap
niAn.winP   = niAn.win*niAn.sRate;
niAn.tStepP = niAn.winP*(1-niAn.pOV);
niAn.winSts = 1:niAn.tStepP:(niAn.numSamp-niAn.winP);
niAn.numWin = length(niAn.winSts);

niAn.sRateDN  = sRate/niAn.dnSamp;
niAn.time     = (0:1/sRate:(r-1)/sRate)';
niAn.pertSig  = squeeze(DAQin(:,1,:));
niAn.sensorFC = squeeze(DAQin(:,2,:));
niAn.sensorFN = squeeze(DAQin(:,3,:));
niAn.sensorP  = squeeze(DAQin(:,4,:));
niAn.audioM   = squeeze(DAQin(:,5,:));
niAn.audioH   = squeeze(DAQin(:,6,:));
niAn.sensorO  = squeeze(DAQin(:,7,:));

niAn.time_audio = dnSampleSmoothSignal(niAn.time, niAn.winP, niAn.numWin, niAn.winSts);
niAn.audioMf0 = signalFrequencyAnalysis(niAn.audioM, niAn.sRate, niAn.numTrial, niAn.winP, niAn.numWin, niAn.winSts);
niAn.audioHf0 = signalFrequencyAnalysis(niAn.audioH, niAn.sRate, niAn.numTrial, niAn.winP, niAn.numWin, niAn.winSts);
niAn.f0b      = (mean(mean(niAn.audioMf0(1:45,:),1)));

niAn.audioMf0_norm = 1200*log2(niAn.audioMf0/niAn.f0b); %Normalize the f0 mic and convert to cents
niAn.audioHf0_norm = 1200*log2(niAn.audioHf0/niAn.f0b); %Normalize the f0 head and convert to cents

niAn.aLimits = [0 4 -100 100];

[B,A] = butter(4, 10/(sRate/2)); %Low-pass filter under 40
niAn.sensorFC = filter(B,A,abs(niAn.sensorFC));
niAn.sensorFN = filter(B,A,abs(niAn.sensorFN));

niAn.time_DN     = dnSampleSignal(niAn.time, niAn.dnSamp); %downSampled Time
niAn.pertSig_DN  = dnSampleSignal(niAn.pertSig, niAn.dnSamp); %downSampled Perturbation Signal
niAn.sensorP_DN  = dnSampleSignal(niAn.sensorP, niAn.dnSamp);
niAn.sensorFC_DN = dnSampleSignal(niAn.sensorFC, niAn.dnSamp);
niAn.sensorFN_DN = dnSampleSignal(niAn.sensorFN, niAn.dnSamp);

[niAn.pertTrig, niAn.idxPert] = findPertTrigs(niAn.time_DN, niAn.pertSig_DN, niAn.sRateDN);
[niAn.presTrig, niAn.idxPres] = findPertTrigs(niAn.time_DN, niAn.sensorP_DN, niAn.sRateDN);
[niAn.fSCTrig, niAn.idxFC]    = findPertTrigs(niAn.time_DN, niAn.sensorFC_DN, niAn.sRateDN);  
[niAn.fSNTrig, niAn.idxFN]    = findPertTrigs(niAn.time_DN, niAn.sensorFN_DN, niAn.sRateDN); 

[niAn.lagsPres, niAn.meanLagTimeP] = calcMeanLags(niAn.pertTrig, niAn.presTrig);
[niAn.lagsFC, niAn.meanLagTimeFC]  = calcMeanLags(niAn.pertTrig, niAn.fSCTrig);
[niAn.lagsFN, niAn.meanLagTimeFN]  = calcMeanLags(niAn.pertTrig, niAn.fSNTrig);

niAn.rangePressures = [];
niAn.timePressures  = [];
for ii = 1:niAn.numTrial
    [dxSmooth, endRiseInd, startFallInd] = findCrossings(niAn.sensorP_DN(:,ii), niAn.tStepP);
    
%     [maxP, maxInd] = max(niAn.sensorP_DN(:,ii));
%     [minP  ] = niAn.sensorP_DN(niAn.idxPert(ii,2), ii);

    onsetPressure  = round(100*niAn.sensorP_DN(endRiseInd,1))/100;
    offsetPressure = round(100*niAn.sensorP_DN(startFallInd,1))/100;
    onsetTime = round(100*niAn.time_DN(endRiseInd))/100;
    offsetTime = round(100*niAn.time_DN(startFallInd))/100;
    
    niAn.timePressures  = cat(1, niAn.timePressures, [onsetTime offsetTime]);
    niAn.rangePressures = cat(1, niAn.rangePressures, [onsetPressure offsetPressure]);
end
niAn.meanRangePressure = mean(niAn.rangePressures, 1);

niAn.riseTimeP = niAn.timePressures(:,1) - niAn.presTrig(:,1);
niAn.meanRiseTimeP = mean(niAn.riseTimeP);
niAn.pLimits = [0 4 0 5];
niAn.fLimits = [0 4 1 5];

niAn.sensorP_Al = alignSensorData(niAn.sRateDN , niAn.numTrial, niAn.idxPert, niAn.sensorP_DN);
niAn.time_Al    = 0:1/niAn.sRateDN :(length(niAn.sensorP_Al)-1)/niAn.sRateDN;
niAn.pLimits_Al = [0 3.5 0 5];


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

function sensorDN = dnSampleSmoothSignal(sensor, winP, numWin, winSts)

sensorDN = [];
for iSt = 1:numWin
    winIdx = winSts(iSt):winSts(iSt) + winP - 1;
    sensorDN = cat(1, sensorDN, mean(sensor(winIdx, :)));
end
end

function sensorf0 = signalFrequencyAnalysis(sensor, fs, numTrial, winP, numWin, winSts)

[B,A] = butter(4, 80/(fs/2), 'high'); %High-pass filter over 60
sensorHP = filter(B,A,abs(sensor), [], 1);

sensorf0 = zeros(numWin, numTrial);
for i = 1:numWin
    for j = 1:numTrial
        winIdx = winSts(i):winSts(i)+ winP - 1;
        sensorf0(i,j) = (calcf0(sensorHP(winIdx, j), fs));
    end
end
end

function f0 = calcf0(x,fs)
% Created by Gabriel Galindo
% Formatted by Dante Smith -12/11/15

lim_inf = ceil(fs/(500));
lim_sup = floor(fs/(50));
U = xcov(x,'unbias');
U = U(ceil(end/2):end);
U = (U(lim_inf:lim_sup)-min(U(lim_inf:lim_sup)))/(max(U(lim_inf:lim_sup)) - min(U(lim_inf:lim_sup)));
[M,P] = findpeaks(U);

if isempty(P)
    f0 = NaN;
else
    P = P(find(M >= 0.9,1,'first'));
    if isempty(P)
        f0 = NaN;
    else
        f0 = fs/(P + lim_inf);
    end

    NFFT = pow2(nextpow2(length(x)/4));
    [Pxx,Fxx] = pwelch(x,NFFT,[],[],fs,'onesided');

    if ~isnan(f0)
        H = Pxx(find(Fxx>=f0,1,'first'));
        if (10*log10(max(Pxx)/H) > 80)
            f0 = NaN;
        end
    end   
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

function [dxSmooth, endRiseInd, startFallInd] = findCrossings(sensor, tStep)

dx = zeros(size(sensor));
for gg = 2:length(sensor)
    dx(gg) = 10000*(sensor(gg) - sensor(gg-1))/tStep;
end
dxSmooth = smooth(dx);
incInds = find(dxSmooth > 0.1);
decInds = find(dxSmooth < -0.2);

incNonIncre = find((diff(incInds) == 1) == 0); %The indices where the values stop increasing
endRiseInd = incInds(incNonIncre(1));

decNonIncre = find((diff(decInds) == 1) == 1); %The indices where the values start decreasing
startFallInd = decInds(decNonIncre(1));

end