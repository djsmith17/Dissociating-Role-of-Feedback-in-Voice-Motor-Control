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

collectNewData         = 1; %Boolean
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
dirs.savFileDir = fullfile(dirs.savFileDir, [expParam.subject '_SensorData.mat']);

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
    
    ForceSensorData.expParam    = expParam;
    ForceSensorData.dirs        = dirs;
    ForceSensorData.DAQin       = DAQin;

    save(dirs.savFileDir, 'ForceSensorData')
else
    load(dirs.savFileDir)
end

pts = length(DAQin);
time = 0:1/expParam.sRateQ:(pts-1)/expParam.sRateQ;
trigs = findPertTrigs(time, DAQin(:,1,:));

[B,A] = butter(4, 40/(expParam.sRateQ/2)); %Low-pass filter under 40

fSensorC = squeeze(DAQin(:,2,:));
fSensorN = squeeze(DAQin(:,3,:));
pSensor  = squeeze(DAQin(:,4,:));

fSensorC  = filter(B,A,abs(fSensorC));
fSensorN  = filter(B,A,abs(fSensorN));

pLimits = [0 4 0 4];
fLimits = [0 4 1 5];
drawDAQsignal(time, fSensorC, fSensorN, pSensor, trigs, pLimits, fLimits, expParam.subject, dirs.savResultsDir, sv2F)
end

function niAn = nidaqAnalysis(expParam, DAQin)

pts = length(DAQin);
time = 0:1/expParam.sRateQ:(pts-1)/expParam.sRateQ;

pert     = squeeze(DAQin(:,1,:));
fSensorC = squeeze(DAQin(:,2,:));
fSensorN = squeeze(DAQin(:,3,:));
pSensor  = squeeze(DAQin(:,4,:));

[B,A] = butter(4, 40/(expParam.sRateQ/2)); %Low-pass filter under 40

niAn.time = time;

end

function trigs = findPertTrigs(time, pertCh)
pertCh = round(squeeze(pertCh)); %Should be step function 0V or 3V
[~, c] = size(pertCh);

trigs = [];
for i = 1:c
    I = find(pertCh(:,i) > 0);
    trigSt = time(I(1));
    trigSp = time(I(end));

    trigs = cat(1, trigs, [trigSt trigSp]);
end
end

function lags = lagCalc(trigs)

end