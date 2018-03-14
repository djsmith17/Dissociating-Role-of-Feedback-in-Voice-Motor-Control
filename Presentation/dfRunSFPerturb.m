function dfRunSFPerturb()
% dfRunSFPerturb()
% Laryngeal Perturbation experiment. This script records acoustic output 
% from a participant as they have their larynx physically displaced.
% NIDAQ signal provides Perturbatron stimulus and Audapter collects and
% manages the recorded acoustic data.
%
% This script calls the following 8 functions:
% dfDirs.m
% dfInitNIDAQ.m
% dfSetAudFB.m
% dfSetTrialOrder.m
% dfMakePertSignal.m
% dfSetVisFB.m
% dfSaveRawData.m
% dfCalcMeanRMS.m
% dfUpdateVisFB.m
%
% This uses the toolbox from MATLAB-Toolboxes
% speechres
%
% This script is also dependent on the following Mathworks Toolboxes
% Signal-Processing Toolbox

close all;
ET = tic;
rng('shuffle');

% Main Experimental prompt: Subject/Run Information
prompt = {'Subject ID:',...
          'Session ID:',...
          'Baseline Loudness (dB SPL):',...
          'Gender ("male" or "female"):',...
          'Balloon:', ...
          'Tightness (inches):'};
name = 'Subject Information';
numlines = 1;
defaultanswer = {'null', 'SF1', '60', 'female', '2.0E_1', 'N/A'};
ExpPrompt = inputdlg(prompt, name, numlines, defaultanswer);

if isempty(ExpPrompt)
    return
end

% Dialogue box asking for what type of Auditory Feedback
AudFB = questdlg('What type of Auditory Feedback?','Auditory Feedback', 'Voice Not Shifted', 'Voice Shifted', 'Masking Noise', 'Masking Noise');
switch AudFB
    case 'Voice Not Shifted'
        AudFBSw = 0;
    case 'Voice Shifted'
        AudFBSw = 1;
    case 'Masking Noise'
        AudFBSw = 2;
end

% Dialogue box asking if Practice set or Full set of trials
num_trials = questdlg('Practice or Full?','Length','Practice','Full','Full');
switch num_trials
    case 'Practice'
        numTrials = 4;
        perCatch  = 1.00;
    case 'Full'
        numTrials = 10;
        perCatch  = 0.50;
end

%Experiment Configurations
expParam.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType      = 'Somatosensory Perturbation_Perceptual';
expParam.subject      = ExpPrompt{1};
expParam.run          = ExpPrompt{2};
expParam.curSess      = [expParam.subject expParam.run];
expParam.targRMS      = str2double(ExpPrompt{3});
expParam.gender       = ExpPrompt{4};
expParam.balloon      = ExpPrompt{5};
expParam.tightness    = ExpPrompt{6};
expParam.niDev        = 'Dev2';              % NIDAQ Device Name. For more information, see dfInitNIDAQ
expParam.trialLen     = 4;                   % Seconds
expParam.numTrial     = numTrials;
expParam.curTrial     = [];
expParam.perCatch     = perCatch;
expParam.AudFB        = AudFB;
expParam.AudFBSw      = AudFBSw;
expParam.AudPert      = '-100 cents ramped'; % Var not used here. Just saving for balance
expParam.AudPertSw    = 1;                   % Var not used here. Just saving for balance
expParam.bVis         = 0;

%Set our dirs based on the project
dirs = dfDirs(expParam.project);
% Folder paths to save data files
dirs.RecFileDir  = fullfile(dirs.RecData, expParam.subject, expParam.run);
dirs.RecWaveDir  = fullfile(dirs.RecFileDir, 'wavFiles');

if exist(dirs.RecFileDir, 'dir') == 0
    mkdir(dirs.RecFileDir)
end
if exist(dirs.RecWaveDir, 'dir') == 0
    mkdir(dirs.RecWaveDir)
end

%Paradigm Configurations
expParam.sRate              = 48000;  % Hardware sampling rate (before downsampling)
expParam.frameLen           = 96;     % Before downsampling
expParam.downFact           = 3;
expParam.sRateAnal          = expParam.sRate/expParam.downFact; %Everything get automatically downsampled! So annoying
expParam.frameLenDown       = expParam.frameLen/expParam.downFact;
expParam.audioInterfaceName = 'MOTU MicroBook'; %'ASIO4ALL' 'Komplete'

%Set up Audapter
Audapter('deviceName', expParam.audioInterfaceName);
Audapter('setParam', 'downFact', expParam.downFact, 0);
Audapter('setParam', 'sRate', expParam.sRateAnal, 0);
Audapter('setParam', 'frameLen', expParam.frameLenDown, 0);
p = getAudapterDefaultParams(expParam.gender);

%Set up Parameters to control NIDAQ and Perturbatron
[s, niCh, nVS]  = dfInitNIDAQ(expParam.niDev, expParam.trialLen);
expParam.sRateQ = s.Rate; % NIDAQ sampling rate
expParam.niCh   = niCh;   % Structure of Channel Names

%Set up OST and PCF Files
expParam.ostFN = fullfile(dirs.Prelim, 'SFPerturbOST.ost'); check_file(expParam.ostFN);
expParam.pcfFN = fullfile(dirs.Prelim, 'SFPerturbPCF.pcf'); check_file(expParam.pcfFN);

%Set up Auditory Feedback (Masking Noise, Pitch-Shift?)
[expParam, p]      = dfSetAudFB(expParam, dirs, p);

% Set up the order of trials (Order of perturbed, control, etc)
expParam.trialType = dfSetTrialOrder(expParam.numTrial, expParam.perCatch);

% Select the trigger points for perturbation onset and offset and creating
% the digital signal to be sent to the NIDAQ
[expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType);

expParam.cuePause  = 1.0; % How long the cue period lasts
expParam.resPause  = 2.0; % How long the rest/VisFB lasts
expParam.boundsRMS = 3;  %+/- dB

% This is where the fun begins
fprintf('\nStarting Trials\n\n')

% Dim the lights (Set the visual Feedback)
[anMsr, H1, H2, H3, fbLines, rec, trigCirc] = dfSetVisFB(expParam.targRMS, expParam.boundsRMS);

%Open the curtains
pause(5);                % Let them breathe a sec
set(H3,'Visible','off'); % Turn off 'Ready?'

DAQin = []; rawData = [];
pltStr = [];
presH = initLiveResult(expParam, 1);
for ii = 1:expParam.numTrial
    expParam.curTrialNum  = ii;
    expParam.curTrial     = ['Trial' num2str(ii)];
    expParam.curSessTrial = [expParam.subject expParam.run expParam.curTrial];
    
    %Used later in audio version
    audStimP = [];
        
    %Set the OST and PCF functions
    Audapter('ost', expParam.ostFN, 0);
    Audapter('pcf', expParam.pcfFN, 0);
    
    %Setup which perturb file we want
    NIDAQsig = [expParam.sigs(:,ii) nVS];
    queueOutputData(s, NIDAQsig);
    
    %Cue to begin trial
    set(H1,'Visible','on');
    pause(expParam.cuePause)
    
    %Phonation Start
    set(H1,'Visible','off');
    set([H2 trigCirc],'Visible','on');
    
    fprintf('Trial %d\n',ii)
    AudapterIO('init', p);
    Audapter('reset');
    Audapter('start');
    
    %Play out the Analog Perturbatron Signal. This will hold script for as
    %long as vector lasts. In this case, 4.0 seconds. 
    [dataDAQ, ~] = s.startForeground;
     
    %Phonation End
    Audapter('stop');
    set([H2 trigCirc],'Visible','off');
    
    % Load the Audapter saved data and save some as wav Files
    data    = dfSaveRawData(expParam, dirs);
    DAQin   = cat(3, DAQin, dataDAQ);
    rawData = cat(1, rawData, data);
       
    %Grab smooth RMS trace from 'data' structure
    rmsMean = dfCalcMeanRMS(data);
    %Compare against baseline and updated Visual Feedback
    [color, newPos] = dfUpdateVisFB(anMsr, rmsMean);

    set(rec, 'position', newPos);
    set(rec, 'Color', color); set(rec, 'FaceColor', color);
    set([rec fbLines], 'Visible', 'on');
    
    pltStr = updateLiveResult(dataDAQ, expParam, pltStr);
    pause(expParam.resPause)
    set([rec fbLines], 'Visible', 'off');
end
close all;
elapsed_time = toc(ET)/60;   % Elapsed Time of the experiment
fprintf('\nElapsed Time: %f (min)\n', elapsed_time)

% Store all the variables and data from the session in a large structure
expParam.elapsedTime = elapsed_time;
DRF.dirs        = dirs;
DRF.expParam    = expParam;
DRF.p           = p;
DRF.audStimP    = audStimP;
DRF.DAQin       = DAQin;
DRF.rawData     = rawData; 

% Save the large data structure (only if not practice trials)
dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.subject expParam.run dirs.saveFileSuffix 'DRF.mat']);
switch num_trials
    case 'Practice'
        return
    case 'Full'
        save(dirs.RecFileDir, 'DRF'); %Only save if it was a full set of trials
end

% qRes = dfQuickAnalysisPlot(DRF)

%Draw the OST progression, if you want to
if expParam.bVis == 1
    OST_MULT = 500; %Scale factor for OST
    visSignals(data, 16000, OST_MULT, savedWavdir)
end
end

function presH = initLiveResult(expParam, defMon)

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

presH = figure('NumberTitle', 'off', 'Color', [1 1 1], 'Position', winPos);

mark = plot([1 1], [-1 5], 'k-', 'LineWidth', 2);
axis([0 3.5 -0.5 5.0])
box off
set(gca,'FontSize', 12,...
        'XTickLabel', {'-1.0' '-0.5' '0' '0.5' '1.0' '1.5' '2.0' '2.5'},...
        'FontWeight', 'bold')
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold') 
ylabel('Pressure (psi)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k') 
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

function visSignals(data, fs, OST_MULT, savedResdir)
plotpos = [200 100];
plotdim = [1000 700];
spectComp = figure('Color', [1 1 1]);
set(spectComp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

frameDur = data.params.frameLen / data.params.sr;
tAxis = 0:frameDur:frameDur * (size(data.rms, 1) - 1);

subplot(2,1,1)
show_spectrogram(data.signalIn, fs, 'noFig');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Microphone In')
box off
plot(tAxis, data.ost_stat * OST_MULT, 'k-');
legend({sprintf('OST status * %d', OST_MULT)});

subplot(2,1,2)
show_spectrogram(data.signalOut, fs, 'noFig');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Headphone Out')
box off

suptitle('ONLINE ''AFA'' TIME WARP SPECTRUM')

plTitle = 'Online AFA Time Warp Spectrum';
saveFileName = [savedResdir plTitle '.png'];

% export_fig(saveFileName)
end