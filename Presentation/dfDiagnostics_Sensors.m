function dfDiagnostics_Sensors(varargin)
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

prompt = {'Subject ID:',...
          'Session ID:',...
          'Number of Trials:',...
          'Percent Perturbed (Dec)',...
          'Collect New Data?:'};
name = 'Subject Information';
numlines = 1;
defaultanswer = {'null', 'DS1', '5', '1', 'yes'};
answer = inputdlg(prompt, name, numlines, defaultanswer);

if isempty(answer)
    return
end

expParam.project       = 'NIDAQSensorDiagnostics';
expParam.expType       = 'Somatosensory Perturbation_Perceptual';
expParam.subject       = answer{1}; %Subject#, Pilot#, null
expParam.run           = answer{2};
expParam.curSess       = [expParam.subject ' ' expParam.run];
expParam.numTrial      = str2double(answer{3});
expParam.perCatch      = str2double(answer{4});
expParam.trialLen      = 4; %Seconds
expParam.niDev         = 'Dev2';
sv2F                   = 1; %Boolean
collectNewData         = answer{5};

dirs = dfDirs(expParam.project);

dirs.RecFileDir    = fullfile(dirs.RecData, expParam.subject, expParam.run);
dirs.SavResultsDir = fullfile(dirs.RecData, expParam.subject, expParam.run); %Where to save results 

if exist(dirs.RecFileDir, 'dir') == 0
    mkdir(dirs.RecFileDir)
end
if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end
dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.subject expParam.run 'NSD.mat']);

if strcmp(collectNewData, 'yes')
    expParam.sRate       = 48000;
    expParam.downFact    = 3;
    expParam.sRateAnal   = expParam.sRate/expParam.downFact; %Everything get automatically downsampled! So annoying
    expParam.resPause    = 3;
    
    [s, niCh, nVS]  = initNIDAQ(expParam.niDev, expParam.trialLen);
    expParam.sRateQ = s.Rate;
    expParam.niCh   = niCh;

    expParam.trialType              = dfSetTrialOrder(expParam.numTrial, expParam.perCatch);
    [expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType, 1);  

    DAQin = []; DAQtime = [];
    for ii = 1:expParam.numTrial
        NIDAQsig = [expParam.sigs(:,ii) nVS];
        queueOutputData(s, NIDAQsig);
        fprintf('Running Trial %d\n', ii)
        
        tic
        [data_DAQ, time] = s.startForeground;
        toc

        DAQin   = cat(3, DAQin, data_DAQ);
        DAQtime = cat(3, DAQtime, time);
        
        pause(expParam.resPause)      
    end
    
    NSD.expParam    = expParam;
    NSD.dirs        = dirs;
    NSD.DAQin       = DAQin;

    save(dirs.RecFileDir, 'NSD')
else
    load(dirs.RecFileDir)
end

[niAn, niRes] = dfAnalysisNIDAQ(dirs, NSD.expParam, NSD.DAQin, 0);

% drawDAQsignal(niAn, 1, dirs.SavResultsDir, sv2F)
drawDAQcombined(niRes, dirs.SavResultsDir, sv2F)
drawDAQAll(niAn, dirs.SavResultsDir, sv2F)
drawDAQPresMic(niAn, dirs.SavResultsDir, sv2F)
end