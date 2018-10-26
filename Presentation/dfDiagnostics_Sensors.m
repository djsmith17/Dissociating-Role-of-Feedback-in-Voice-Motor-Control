function dfDiagnostics_Sensors()
% dfDiagnostics_Sensors() collects data from the NIDAQ in the way that the main experimental A quick test of the force sensors before running the actual experiment.
% This makes sure that the sensors are working they should be and we can
% continue with the experiment. Eventually this will also include the
% pressure sensor. 
%
% This script calls the following (4) functions:
% dfDirs.m
% dfInitNIDAQ.m
% dfMakePertSignal.m
% dfAnalysisNIDAQ.m
% drawDAQsignal.m
close all;

% Main Experimental prompt: Subject/Run Information
prompt = {'Subject ID:',...
          'Session ID:',...
          'Number of Trials:',...
          'Percent Perturbed (Dec)',...
          'Balloon:', ...
          'Collect New Data?:'};
name = 'Subject Information';
numlines = 1;
defaultanswer = {'null', 'DS1', '4', '1', '2.0K_4','yes'};
ExpPrompt = inputdlg(prompt, name, numlines, defaultanswer);

if isempty(ExpPrompt)
    return
end

%Experiment Configurations
expParam.project       = 'NIDAQSensorDiagnostics';
expParam.expType       = 'Somatosensory Perturbation_Perceptual';
expParam.subject       = ExpPrompt{1}; %Subject#, Pilot#, null
expParam.run           = ExpPrompt{2};
expParam.curSess       = [expParam.subject ' ' expParam.run];
expParam.gender        = 'N/A';
expParam.niDev         = 'Dev2';                  % NIDAQ Device Name. For more information, see dfInitNIDAQ
expParam.trialLen      = 4;                       % Seconds
expParam.numTrial      = str2double(ExpPrompt{3});
expParam.curTrial      = [];
expParam.perCatch      = str2double(ExpPrompt{4});
expParam.balloon       = ExpPrompt{5};
expParam.AudFB         = 'Masking Noise';
expParam.AudFBSw       = 2;
expParam.resPause      = 6;
expParam.trialLenLong  = expParam.numTrial*(expParam.trialLen + expParam.resPause);
expParam.sigLong       = [];

sv2F                   = 1; %Boolean
collectNewData         = ExpPrompt{6};

%Set our dirs based on the project
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
    expParam.sRateAnal   = expParam.sRate/expParam.downFact;
    
    %Set up Parameters to control NIDAQ and Perturbatron
    [s, niCh, nVS]  = dfInitNIDAQ(expParam.niDev, expParam.trialLen);
    expParam.sRateQ = s.Rate; % NIDAQ sampling rate
    expParam.niCh   = niCh;   % Structure of Channel Names

    % Set up the order of trials (Order of perturbed, control, etc)
    expParam.trialType = dfSetTrialOrder(expParam.numTrial, expParam.perCatch);
    
    % Select the trigger points for perturbation onset and offset and creating
    % the digital signal to be sent to the NIDAQ
    [expParam.sigs, expParam.trigs, expParam.vSigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, 1);  
    
    DAQin = []; DAQtime = [];
    LR = LiveSensorResult(expParam, 2);
    for ii = 1:expParam.numTrial
        expParam.curTrialNum  = ii;
        expParam.curTrial     = ['Trial' num2str(ii)];
        expParam.curSessTrial = [expParam.subject expParam.run expParam.curTrial];
        
        %Setup which perturb file we want
        NIDAQsig = [expParam.sigs(:,ii) expParam.vSigs(:,ii)];
        queueOutputData(s, NIDAQsig);        
        
        fprintf('Running Trial %d\n', ii)
        % Play out the Analog Perturbatron Signal. 
        % This will hold script for as long as vector lasts (4s) 
        [data_DAQ, time] = s.startForeground;

        DAQin   = cat(3, DAQin, data_DAQ);
        DAQtime = cat(3, DAQtime, time);
        
        LR = LR.updateLiveResult(data_DAQ, ii);
        pause(expParam.resPause)      
    end
    
    % Store all the variables and data from the session in a large structure
    NSD.expParam    = expParam;
    NSD.dirs        = dirs;
    NSD.DAQin       = DAQin;

    % Save the large data structure
    save(dirs.RecFileDir, 'NSD')
else
    load(dirs.RecFileDir)
end

% f0b = 100;
% pF  = 1;
% iRF = 0;
% [niAn, niRes] = dfAnalysisNIDAQ(dirs, NSD.expParam, NSD.DAQin, f0b, 0, iRF, pF);
% 
% niRes.numPertTrialsNi = niRes.numPertTrials;
% 
% % drawDAQsignal(niAn, 1, dirs.SavResultsDir, sv2F)
% drawDAQAlignedPressure(niRes, dirs.SavResultsDir, sv2F)
% drawDAQAll(niAn, dirs.SavResultsDir, sv2F)
% drawDAQPresMic(niAn, dirs.SavResultsDir, sv2F)
end