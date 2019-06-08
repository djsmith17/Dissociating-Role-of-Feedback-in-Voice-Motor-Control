function dfDiagnostics_Sensors()
% dfDiagnostics_Sensors() collects data from the NIDAQ under standard 
% conditions to test how the recordings look before entering a full
% experiment.
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
defaultanswer = {'Test', 'DS1', '3', '1', 'DefaultBalloon','yes'};
ExpPrompt = inputdlg(prompt, name, numlines, defaultanswer);

sensorPType = 'Seven';
doubleSens  = 0;

if isempty(ExpPrompt)
    returns
end

%Experiment Configurations
expParam.project       = 'NIDAQSensorDiagnostics';
expParam.expType       = 'Somatosensory Perturbation_Perceptual';
expParam.subject       = ExpPrompt{1}; %Subject#, Pilot#, null
expParam.run           = ExpPrompt{2};
expParam.curSess       = [expParam.subject ' ' expParam.run];
expParam.gender        = 'N/A';
expParam.niDev         = 'Dev1';                  % NIDAQ Device Name. For more information, see dfInitNIDAQ
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
expParam.sensorPType   = sensorPType;

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
    [expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, 1);  
    
    DAQin = []; DAQtime = [];
    LR = LiveSensorResult(expParam, 2);
    for ii = 1:expParam.numTrial
        expParam.curTrialNum  = ii;
        expParam.curTrial     = ['Trial' num2str(ii)];
        expParam.curSessTrial = [expParam.subject expParam.run expParam.curTrial];
        
        %Setup which perturb file we want
        NIDAQsig = [expParam.sigs(:,ii) nVS];
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

f0b = 100;
pF  = 1;
aD = 0;
[niAn, niRes] = dfAnalysisNIDAQ(dirs, NSD.expParam, NSD.DAQin, f0b, 0, aD, pF);

niRes.numPertTrialsNi = niRes.numPertTrials;

% drawDAQsignal(niAn, 1, dirs.SavResultsDir, sv2F)
drawDAQAlignedPressure(niRes, dirs.SavResultsDir, sv2F, 0)
drawRawVoltagePressure(niAn, dirs.SavResultsDir, 0)
if doubleSens == 1
    drawDAQAlignedPressure(niRes, dirs.SavResultsDir, sv2F, 1)
    drawRawVoltagePressure(niAn, dirs.SavResultsDir, 1)
    CompareVoltagesInternalExternal()
end
% drawDAQAll(niAn, dirs.SavResultsDir, sv2F)
% drawDAQPresMic(niAn, dirs.SavResultsDir, sv2F)
end

function drawRawVoltagePressure(niAn, saveResultsDir, intFlag)
curSess = niAn.curSess;
plotpos = [500 300];
plotdim = [850 600];
rawV = figure('Color', [1 1 1]);
set(rawV, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

if intFlag == 1
    sensor = niAn.sensorFN;
    suffix = '_Int';
else
    sensor = niAn.sensorP;
    suffix = '';
end

trialColors = distinguishable_colors(niAn.numTrial);
maxVals = [];
for i = 1:niAn.numTrial
    plot(niAn.time, sensor(:,i), 'Color', trialColors(i,:))
    maxVal = max(sensor(:,i));
    maxVals = cat(1, maxVals, maxVal);
    hold on
end
xlabel('Time(s)','FontSize', 12, 'FontWeight', 'bold')
ylabel('Voltage (V)','FontSize', 12, 'FontWeight', 'bold')
title({curSess; 'Raw Voltage'},'FontSize', 12, 'FontWeight', 'bold')
axis([0 4 0 5])
box off

meanMaxVal = round(mean(maxVals), 2);
annotation('TextBox', [0.65 0.8 0.1 0.1],...
           'String', ['Mean Plateau Val: ' num2str(meanMaxVal) 'V'],...
           'LineStyle','none',...
            'FontWeight','bold',...
            'FontSize',15,...
            'FontName','Arial')
        
plTitle = [curSess  '_RawVoltageTrace' suffix '.jpg'];     
saveFileName = fullfile(saveResultsDir, plTitle);
export_fig(saveFileName) 
end