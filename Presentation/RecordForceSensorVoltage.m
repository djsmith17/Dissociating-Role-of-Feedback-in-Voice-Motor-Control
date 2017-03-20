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
close all;

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

dirs.savFileDir    = fullfile(dirs.RecData, expParam.subject);
dirs.savResultsDir = fullfile(dirs.RecData, expParam.subject); %Where to save results 

if exist(dirs.savFileDir, 'dir') == 0
    mkdir(dirs.savFileDir)
end
if exist(dirs.savResultsDir, 'dir') == 0
    mkdir(dirs.savFileDir)
end

expParam.sRate       = 48000;
expParam.downFact    = 3;
expParam.sRateAnal   = expParam.sRate/expParam.downFact; %Everything get automatically downsampled! So annoying

[s, niCh, nVS] = initNIDAQ(expParam.trialLen, 'Dev3');
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

limits  = [0 4 0 5];
drawDAQsignal(expParam.sRateQ, expParam.trigs(:,:,1), DAQin, limits, expParam.subject, dirs.savResultsDir)

ForceSensorData.expParam    = expParam;
ForceSensorData.dirs        = dirs;
ForceSensorData.DAQin       = DAQin;

dirs.savFileDir = fullfile(dirs.savFileDir, [expParam.subject '_SensorData.mat']);
save(dirs.savFileDir, 'ForceSensorData')
end