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
defaultanswer = {'null', 'DS1', '2', '1', '2.0K_4','yes'};
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
expParam.resPause      = 3;
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
    [expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, 1);  
    
    DAQin = []; DAQtime = [];
    pltStr = [];
    presH = initLiveResult(expParam, 2);
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
        
        pltStr = updateLiveResult(data_DAQ, expParam, pltStr);
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
iRF = 0;
[niAn, niRes] = dfAnalysisNIDAQ(dirs, NSD.expParam, NSD.DAQin, f0b, 0, iRF, pF);

niRes.numPertTrialsNi = niRes.numPertTrials;

% drawDAQsignal(niAn, 1, dirs.SavResultsDir, sv2F)
drawDAQAlignedPressure(niRes, dirs.SavResultsDir, sv2F)
% drawDAQAll(niAn, dirs.SavResultsDir, sv2F)
% drawDAQPresMic(niAn, dirs.SavResultsDir, sv2F)
end

function resH = initLiveResult(expParam, defMon)
% resH = initLiveResult(expParam, defMon) creates a figure by which to
% display recorded signal during a recording session. This is currently
% configured for recording pressure in the perturbatron balloon.
% 
% expParam: Experimental paramenters of the recording
% defMon  : Monitor definitions
%
% resH    : Figure handle for the generated figure. 

curSess  = expParam.curSess;
balloon  = expParam.balloon;
balloon(strfind(balloon, '_')) = '';

monitorSize = get(0, 'Monitor');
numMon = size(monitorSize, 1);
plotDim = [800 600];

if numMon == 2 && defMon == 2
    [~, mon] = max(monitorSize(:,1));
    
    halfW  = monitorSize(mon, 3)/2;
    halfWD = halfW - plotDim(1)/2 + monitorSize(mon, 1) - 1;
    
    figPosition = [halfWD 80 plotDim];
else
    
    halfW = monitorSize(1, 3)/2;
    halfWD = halfW - plotDim(1)/2 + monitorSize(1, 1) - 1;
    
    figPosition = [halfWD 80 plotDim];
end
winPos = figPosition;

resH = figure('NumberTitle', 'off', 'Color', [1 1 1], 'Position', winPos);

mark = plot([1 1], [-1 5], 'k-', 'LineWidth', 2);
axis([0 3.5 -0.5 5.0])
box off
set(gca,'FontSize', 12,...
        'XTickLabel', {'-1.0' '-0.5' '0' '0.5' '1.0' '1.5' '2.0' '2.5'},...
        'FontWeight', 'bold')
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold') 
ylabel('Voltage (V)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k') 
title({'Pressure Recording, Live Result';
       curSess;
       ['Balloon: ' balloon]})

hold on
end

function pltStr = updateLiveResult(daqIn, expParam, pltStr)

sig      = daqIn(:,4);
numTrial = expParam.numTrial;
curTrial = expParam.curTrialNum;
fs       = expParam.sRateQ;
trigs    = expParam.trigs(:,:,2);
trialColors = distinguishable_colors(numTrial);

St = trigs(curTrial,1) - fs*1 + 1;
Sp = trigs(curTrial,1) + fs*2.5;

sigSnip = sig(St:Sp);
time    = (0:1/fs :(length(sigSnip)-1)/fs)';

tag = ['Trial ' num2str(curTrial)];
trPrs = plot(time, sigSnip, 'LineWidth', 2, 'Color', trialColors(curTrial, :));

if curTrial == 1
    pltStr.tag = {tag};
    pltStr.curve = trPrs;
else
    pltStr.tag   = cat(1, pltStr.tag, tag);
    pltStr.curve = cat(1, pltStr.curve, trPrs);
end

lgd = legend(pltStr.curve, pltStr.tag);
set(lgd, 'box', 'off',...
         'location', 'NorthWest'); 
end