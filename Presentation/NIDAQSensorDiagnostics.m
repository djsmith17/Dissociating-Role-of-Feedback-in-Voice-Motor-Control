function NIDAQSensorDiagnostics(varargin)
%A quick test of the force sensors before running the actual experiment.
%This makes sure that the sensors are working they should be and we can
%continue with the experiment. Eventually this will also include the
%pressure sensor. 

%This script calls the following (4) functions:
%sfDirs.m
%initNIDAQ.m
%createPerturbSignal.m
%drawDAQsignal.m
close all;

if isempty(varargin)
    numTrial = 4; 
else
    numTrial = varargin{1};
end

collectNewData         = 0; %Boolean
sv2F                   = 1; %Boolean

expParam.project       = 'NIDAQSensorDiagnostics';
expParam.expType       = 'Somatosensory Perturbation_Perceptual';
expParam.subject       = 'Pilot0_noPhonation'; %Subject#, Pilot#, null
expParam.numTrial      = numTrial; %Experimental trials = 40
expParam.trialLen      = 4; %Seconds

dirs = sfDirs(expParam.project);

dirs.savFileDir    = fullfile(dirs.RecData, expParam.subject);
dirs.savResultsDir = fullfile(dirs.RecData, expParam.subject); %Where to save results 

if exist(dirs.savFileDir, 'dir') == 0
    mkdir(dirs.savFileDir)
end
if exist(dirs.savResultsDir, 'dir') == 0
    mkdir(dirs.savFileDir)
end
dirs.savFileDir = fullfile(dirs.savFileDir, [expParam.subject '_NSD.mat']);

if collectNewData == 1
    expParam.sRate       = 48000;
    expParam.downFact    = 3;
    expParam.sRateAnal   = expParam.sRate/expParam.downFact; %Everything get automatically downsampled! So annoying

    [s, niCh, nVS]  = initNIDAQ(expParam.trialLen, 'Dev3');
    expParam.sRateQ = s.Rate;
    expParam.niCh   = niCh;

    expParam.trialType = ones(expParam.numTrial,1);

    [expParam.sigs, expParam.trigs] = createPerturbSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType);

    expParam.resPause = 1;

    DAQin = [];
    for ii = 1:expParam.numTrial
        NIDAQsig = [expParam.sigs(:,ii) nVS];
        queueOutputData(s, NIDAQsig);
        fprintf('Running Trial %d\n', ii)
        [data_DAQ, time] = s.startForeground;

        DAQin = cat(3, DAQin, data_DAQ);

        pause(expParam.resPause)      
    end
    
    NSD.expParam    = expParam;
    NSD.dirs        = dirs;
    NSD.DAQin       = DAQin;

    save(dirs.savFileDir, 'NSD')
else
    load(dirs.savFileDir)
end

niAn = nidaqAnalysis(NSD.expParam, NSD.DAQin);

pLimits = [0 4 0 4];
fLimits = [0 4 1 5];
drawDAQsignal(niAn.time, niAn.fSensorC, niAn.fSensorN, niAn.pSensor, niAn.trigs, niAn, pLimits, fLimits, NSD.expParam.subject, dirs.savResultsDir, sv2F)

pLimits = [0 3.5 0 4];
drawDAQcombined(niAn.timeAl, niAn.pSensorAl, niAn.trigs, niAn, pLimits, NSD.expParam.subject, dirs.savResultsDir, sv2F)
end

function niAn = nidaqAnalysis(expParam, DAQin)

sRate = expParam.sRateQ;

[r, c, n] = size(DAQin);
time = 0:1/sRate:(r-1)/sRate;

numTrial = n; 
pert     = squeeze(DAQin(:,1,:));
fSensorC = squeeze(DAQin(:,2,:));
fSensorN = squeeze(DAQin(:,3,:));
pSensor  = squeeze(DAQin(:,4,:));

[B,A] = butter(4, 40/(sRate/2)); %Low-pass filter under 40
fSensorC  = filter(B,A,abs(fSensorC));
fSensorN  = filter(B,A,abs(fSensorN));

[pertTrig, pertThresh, pertidx] = findPertTrigs(time, pert);
[presTrig, presThresh, presidx] = findPertTrigs(time, pSensor);
[fSCTrig, fSCThresh, fSCidx]   = findPertTrigs(time, fSensorC);  
[fSNTrig, fSNThresh, fSNidx]   = findPertTrigs(time, fSensorN); 

[PresLags, PresLagVals] = calcMeanLags(pertTrig, presTrig);
[fSCLags, fSCLagVals]   = calcMeanLags(pertTrig, fSCTrig);
[fSNLags, fSNLagVals]   = calcMeanLags(pertTrig, fSNTrig);

rangePressures = [];
for ii = 1:numTrial
    onsetPressure  = round(100*max(pSensor(:,ii)))/100;
    offsetPressure = round(100*pSensor(pertidx(ii,2), ii))/100;
    rangePressures = cat(1, rangePressures, [onsetPressure offsetPressure]);
end

pSensorAl = alignSensorData(numTrial, sRate, pSensor, pertidx);
timeAl = 0:1/sRate:(length(pSensorAl)-1)/sRate;

niAn.numTrial = numTrial;
niAn.time = time;
niAn.pert = pert;
niAn.fSensorC = fSensorC;
niAn.fSensorN = fSensorN;
niAn.pSensor  = pSensor;

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

niAn.pSensorAl = pSensorAl;
niAn.timeAl = timeAl;
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

function sensorAl = alignSensorData(numTrial, sRate, sensor, idx)

sensorAl = [];
for ii = 1:numTrial
    St = idx(ii,1) - sRate*1;
    Sp = idx(ii,1) + sRate*2.5;
    
    sensorAl = cat(2, sensorAl, sensor(St:Sp,ii));
end
end