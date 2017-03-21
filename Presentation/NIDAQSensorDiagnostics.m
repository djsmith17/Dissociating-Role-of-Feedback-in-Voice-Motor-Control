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
expParam.subject       = 'null'; %Subject#, Pilot#, null
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
drawDAQsignal(niAn.time, niAn.fSensorC, niAn.fSensorN, niAn.pSensor, niAn.trigs, pLimits, fLimits, NSD.expParam.subject, dirs.savResultsDir, sv2F)
end

function niAn = nidaqAnalysis(expParam, DAQin)

sRate = expParam.sRateQ;

pts = length(DAQin);
time = 0:1/sRate:(pts-1)/sRate;

pert     = squeeze(DAQin(:,1,:));
fSensorC = squeeze(DAQin(:,2,:));
fSensorN = squeeze(DAQin(:,3,:));
pSensor  = squeeze(DAQin(:,4,:));

[B,A] = butter(4, 40/(sRate/2)); %Low-pass filter under 40
fSensorC  = filter(B,A,abs(fSensorC));
fSensorN  = filter(B,A,abs(fSensorN));

[pertTrig, pertThresh] = findPertTrigs(time, pert);
[presTrig, presThresh] = findPertTrigs(time, pSensor);
[fSCTrig, fSCThresh]   = findPertTrigs(time, fSensorC);  
[fSNTrig, fSNThresh]   = findPertTrigs(time, fSensorN); 

[Preslags, PresLagVals] = calcMeanLags(pertTrig, presTrig);
[fSClags, fSCLagVals]   = calcMeanLags(pertTrig, fSCTrig);
[fSNlags, fSNLagVals]   = calcMeanLags(pertTrig, fSNTrig);

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

niAn.PresLags = Preslags;
niAn.PresLagVals = PresLagVals;
niAn.fSCLags  = fSClags;
niAn.fSCLagVals = fSCLagVals;
niAn.fSNLags = fSNlags;
niAn.fSNLagVals = fSNLagVals;
end

function [trigs, threshes] = findPertTrigs(time, pertCh)
pertCh = round(pertCh); %Should be step function 0V or 3V
[~, c] = size(pertCh);

trigs = [];
threshes = [];
for i = 1:c
    thresh = mean(pertCh(2000:4000, i));
    
    I = find(pertCh(:,i) > thresh);
    trigSt = time(I(1));
    trigSp = time(I(end));

    trigs = cat(1, trigs, [trigSt trigSp]);
    threshes = cat(1, threshes, thresh);
end
end

function [lags, lagVals] = calcMeanLags(pertTrig, sensorTrig)

lags = sensorTrig - pertTrig;
lagsMean = mean(lags, 1);
lagsSTD  = std(lags, 0, 2);

SEM = lagsSTD/sqrt(length(lags)); 
CIM = 1.96*SEM;

lagVals = [lagsMean, CIM];
end