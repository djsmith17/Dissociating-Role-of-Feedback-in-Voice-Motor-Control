function dfDiagnostics_SensorsLong(varargin)
%A quick test of the force sensors before running the actual experiment.
%This makes sure that the sensors are working they should be and we can
%continue with the experiment. Eventually this will also include the
%pressure sensor. 

%This script calls the following (4) functions:
%dfDirs.m
%initNIDAQ.m
%dfMakePertSignal.m
%drawDAQsignal.m
close all;

if isempty(varargin)
    numTrial = 10; 
else
    numTrial = varargin{1};
end

collectNewData         = 0; %Boolean
sv2F                   = 1; %Boolean

expParam.project       = 'NIDAQSensorDiagnostics';
expParam.expType       = 'Somatosensory Perturbation_Perceptual';
expParam.subject       = 'BalloonD_EmptyAirLong'; %Subject#, Pilot#, null
expParam.run           = 'Run1';
expParam.curSess       = [expParam.subject ' ' expParam.run];
expParam.numTrial      = numTrial; %Experimental trials = 40
expParam.trialLen      = 4; %Seconds
expParam.perCatch      = 1;

dirs = dfDirs(expParam.project);

dirs.RecFileDir    = fullfile(dirs.RecData, expParam.subject);
dirs.SavResultsDir = fullfile(dirs.RecData, expParam.subject); %Where to save results 

if exist(dirs.RecFileDir, 'dir') == 0
    mkdir(dirs.RecFileDir)
end
if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end
dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.subject '_DiagSensors.mat']);

if collectNewData == 1
    expParam.sRate       = 48000;
    expParam.downFact    = 3;
    expParam.sRateAnal   = expParam.sRate/expParam.downFact; %Everything get automatically downsampled! So annoying
    expParam.resPause    = 1;
    
    %New!
    expParam.trialLenLong = expParam.numTrial*(expParam.trialLen + expParam.resPause);
    
    [s, niCh, nVS]  = initNIDAQ(expParam.trialLenLong, 'Dev2');
    expParam.sRateQ = s.Rate;
    expParam.niCh   = niCh;

    expParam.trialType = dfSetTrialOrder(expParam.numTrial, expParam.perCatch);

    [expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType);

    %New!
    expParam.resPauseVec = zeros(expParam.sRateQ*expParam.resPause,1);
    expParam.sigLong     = [];
    for w = 1:expParam.numTrial
        expParam.sigLong = [expParam.sigLong; expParam.sigs(:,w); expParam.resPauseVec];
    end

    NIDAQsig = [expParam.sigLong, nVS];
    queueOutputData(s, NIDAQsig);
    fprintf('Running Trial %d\n', 1)
    [data_DAQ, time] = s.startForeground;
    
    NSD.expParam    = expParam;
    NSD.dirs        = dirs;
    NSD.DAQin       = data_DAQ;

    save(dirs.RecFileDir, 'NSD')
else
    load(dirs.RecFileDir)
end

niAn = dfAnalysisNIDAQ(NSD.expParam, NSD.DAQin);

niAn.pLimits = [0 50 0 5];
niAn.fLimits = [0 50 1 5];
drawDAQsignal_Long(niAn, 1, dirs.SavResultsDir, sv2F)
end