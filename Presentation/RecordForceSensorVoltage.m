function RecordForceSensorVoltage(varargin)
%A quick test of the force sensors before running the actual experiment.
%This makes sure that the sensors are working they should be and we can
%continue with the experiment. Eventually this will also include the
%pressure sensor. 

%This script calls the following (4) functions:
%sfDirs.m
%initNIDAQ.m
%createPerturbSignal.m
%drawDAQsignal.m

if isempty(varargin)
    numTrial = 4; 
else
    numTrial = varargin{1};
end

expParam.project       = 'Calibration_Force Sensor';
expParam.expType       = 'Somatosensory Perturbation_Perceptual';
expParam.subject       = 'null'; %Subject#, Pilot#, null
expParam.numTrial      = numTrial; %Experimental trials = 40
expParam.trialLen      = 4; %Seconds

dirs = sfDirs(expParam.project);

dirs.saveFileDir    = fullfile(dirs.Data, expParam.subject);
dirs.saveResultsDir = fullfile(dirs.Data, expParam.subject); %Where to save results 

if exist(dirs.saveFileDir, 'dir') == 0
    mkdir(dirs.saveFileDir)
end
if exist(dirs.saveResultsDir, 'dir') == 0
    mkdir(dirs.saveFileDir)
end

expParam.sRate       = 48000;
expParam.downFact    = 3;
expParam.sRateAnal   = expParam.sRate/expParam.downFact; %Everything get automatically downsampled! So annoying

s = initNIDAQ;
expParam.sRateQ = s.Rate;

expParam.trialType = ones(expParam.numTrial,1);

[expParam.sigs, expParam.trigs] = createPerturbSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType);

negVolSrc = zeros(expParam.sRateQ*expParam.trialLen, 1) - 1;
negVolSrc(1) = 0; negVolSrc(end) = 0;

expParam.resPause = 1;

DAQin = [];
for ii = 1:expParam.numTrial
    NIDAQsig = [expParam.sigs(:,ii) negVolSrc];
    queueOutputData(s, NIDAQsig);
    fprintf('Running Trial %d\n', ii)
    [data_DAQ, time] = s.startForeground;
    
    DAQin = cat(3, DAQin, data_DAQ);
    
    pause(expParam.resPause)      
end

drawDAQsignal(expParam.sRateQ, expParam.trigs(:,:,2), DAQin, expParam.subject, dirs.saveResultsDir)

ForceSensorData.expParam    = expParam;
ForceSensorData.dirs        = dirs;
ForceSensorData.DAQin       = DAQin;

dirs.saveFileDir = fullfile(dirs.saveFileDir, [expParam.subject '_SensorData.mat']);
save(dirs.saveFileDir, 'ForceSensorData')
end