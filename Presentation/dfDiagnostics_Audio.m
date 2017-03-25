function dfDiagnostics_Audio(varargin)
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
    numTrial = 4; 
else
    numTrial = varargin{1};
end

collectNewData         = 0; %Boolean
sv2F                   = 1; %Boolean

expParam.project       = 'Diagnostics_Audio';
expParam.expType       = 'Auditory Perturbation_Perceptual';
expParam.subject       = 'null'; %Subject#, Pilot#, null
expParam.numTrial      = numTrial; %Experimental trials = 40
expParam.trialLen      = 4; %Seconds

dirs = dfDirs(expParam.project);

dirs.savFileDir    = fullfile(dirs.RecData, expParam.subject);
dirs.savResultsDir = fullfile(dirs.RecData, expParam.subject); %Where to save results 

if exist(dirs.savFileDir, 'dir') == 0
    mkdir(dirs.savFileDir)
end
if exist(dirs.savResultsDir, 'dir') == 0
    mkdir(dirs.savFileDir)
end
dirs.savFileDir = fullfile(dirs.savFileDir, [expParam.subject '_DiagAud.mat']);

if collectNewData == 1
    expParam.sRate       = 48000;
    expParam.downFact    = 3;
    expParam.sRateAnal   = expParam.sRate/expParam.downFact; %Everything get automatically downsampled! So annoying

    [s, niCh, nVS]  = initNIDAQ(expParam.trialLen, 'Dev3');
    expParam.sRateQ = s.Rate;
    expParam.niCh   = niCh;

    expParam.trialType = ones(expParam.numTrial,1);

    [expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType);

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

niAn = dfAnalysisNIDAQ(NSD.expParam, NSD.DAQin);

pLimits = [0 4 0 4];
fLimits = [0 4 1 5];
drawDAQsignal(niAn.time, niAn.fSensorC, niAn.fSensorN, niAn.pSensor, niAn.trigs, niAn, pLimits, fLimits, NSD.expParam.subject, dirs.savResultsDir, sv2F)

pLimits = [0 3.5 0 4];
drawDAQcombined(niAn.timeAl, niAn.pSensorAl, niAn.trigs, niAn, pLimits, NSD.expParam.subject, dirs.savResultsDir, sv2F)
end