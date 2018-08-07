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
subject    = 'DRF_MN4';   % Subject#, Pilot#, null
run        = prompt4RunName();

balloon    = '2E4';     % Which perturbation balloon?
tightness  = 'n/a';     % (inches of slack in bungie cord)
baseV      = 'BV1';

% Dialogue box asking for what type of Auditory Feedback
AudFB = questdlg('What type of Auditory Feedback?','Auditory Feedback', 'Voice Feedback', 'AC Masking Noise', 'AC/BC Masking Noise', 'AC Masking Noise');
switch AudFB
    case 'Voice Feedback'
        headCk = questdlg('Are the BC Headphones unplugged?','BC Headphones','Yes', 'Yes');
        AudFBSw = 0;
    case 'AC Masking Noise'
        headCk = questdlg('Are the BC Headphones unplugged?','BC Headphones','Yes', 'Yes');
        AudFBSw = 2;
    case 'AC/BC Masking Noise'
        headCk = questdlg('Are the BC Headphones plugged in and on?','BC Headphones','Yes', 'Yes');
        AudFBSw = 2;
end

% Dialogue box asking if Practice set or Full set of trials
recType = questdlg('Practice, Diagnostic, or Full?','Length','Practice', 'Diagnostic', 'Full','Full');
switch recType
    case 'Practice'
        numTrials = 4;
        perCatch  = 1.00;
    case 'Diagnostic'
        numTrials = 10;
        perCatch  = 0.50;
    case 'Full'
        numTrials = 40;
        perCatch  = 0.25;
end

fprintf('\nBeginning %s set of recordings for\n%s %s with %s and Balloon %s\n\n', recType, subject, run, AudFB, balloon)

%Experiment Configurations
expParam.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType      = 'Somatosensory Perturbation_Perceptual';
expParam.subject      = subject;
expParam.run          = run;
expParam.curSess      = [expParam.subject expParam.run];
expParam.balloon      = balloon;
expParam.tightness    = tightness;
expParam.InflaVarNm   = 'N/A';
expParam.niDev        = 'Dev2';              % NIDAQ Device Name. For more information, see dfInitNIDAQ
expParam.trialLen     = 4;                   % Seconds
expParam.numTrial     = numTrials;
expParam.curTrial     = [];
expParam.perCatch     = perCatch;
expParam.headGain     = 5;                   % Output gain above the input
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
dirs.BaseFile    = fullfile(dirs.RecData, expParam.subject, baseV, [expParam.subject baseV 'DRF.mat']);

if exist(dirs.RecFileDir, 'dir') == 0
    mkdir(dirs.RecFileDir)
end
if exist(dirs.RecWaveDir, 'dir') == 0
    mkdir(dirs.RecWaveDir)
end

[expParam.f0b,...
 expParam.targRMS,...
 expParam.rmsB,...
 expParam.gender,...
 expParam.age] = loadBaselineVoice(dirs);

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
expParam.buffPause = 0.8; %Give them a moment to start speaking
expParam.endPause  = 0.5;
expParam.resPause  = 2.0; % How long the rest/VisFB lasts
expParam.boundsRMS = 3;   % +/- dB

% This is where the fun begins
fprintf('\nStarting Trials\n\n')

% Dim the lights (Set the visual Feedback)
[anMsr, H1, H2, H3, fbLines, rec, trigCirc] = dfSetVisFB(expParam.curSess, expParam.targRMS, expParam.boundsRMS);

%Open the curtains
pause(5);                % Let them breathe a sec
set(H3,'Visible','off'); % Turn off 'Ready?'

DAQin = []; rawData = [];
loudResults = [];
LR = LiveSensorResult(expParam, 1);
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
    pause(expParam.buffPause)
    
    %Play out the Analog Perturbatron Signal. This will hold script for as
    %long as vector lasts. In this case, 4.0 seconds. 
    [dataDAQ, ~] = s.startForeground;
     
    %Phonation End
    set([H2 trigCirc],'Visible','off');
    pause(expParam.endPause)
    Audapter('stop');
    
    % Load the Audapter saved data and save as wav Files
    data    = AudapterIO('getData'); % This will need to become a try statement again
    DAQin   = cat(3, DAQin, dataDAQ);
    rawData = cat(1, rawData, data);
       
    %Grab smooth RMS trace from 'data' structure
    rmsMean = dfCalcMeanRMS(data, expParam.rmsB);
    %Compare against baseline and updated Visual Feedback
    [color, newPos, loudResult] = dfUpdateVisFB(anMsr, rmsMean);
    loudResults = cat(1, loudResults, loudResult);
    dispLoudnessResult(loudResult)

    set(rec, 'position', newPos);
    set(rec, 'Color', color); set(rec, 'FaceColor', color);
    set([rec fbLines], 'Visible', 'on');
    
    LR = LR.updateLiveResult(dataDAQ, ii);
    
    dfSaveWavRec(data, expParam, dirs);    
    pause(expParam.resPause)
    set([rec fbLines], 'Visible', 'off');
end
close all;
elapsed_time = toc(ET)/60;   % Elapsed Time of the experiment
fprintf('\nElapsed Time: %f (min)\n', elapsed_time)

% Store all the variables and data from the session in a large structure
expParam.elapsedTime = elapsed_time;
expParam.loudResults = loudResults;
DRF.dirs        = dirs;
DRF.expParam    = expParam;
DRF.p           = p;
DRF.audStimP    = audStimP;
DRF.DAQin       = DAQin;
DRF.rawData     = rawData; 

switch recType
    case 'Practice'
        for ii = 1:expParam.numTrial
            dfShowPraatSpect(dirs, expParam.curSess, ['Trial' num2str(ii)])
        end
%         DRF.qRes = dfAnalysisAudioQuick(DRF, 1);
end

% Save the large data structure (only if not practice trials)
dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.subject expParam.run dirs.saveFileSuffix 'DRF.mat']);
fprintf('\nSaving recorded data at:\n%s\n\n', dirs.RecFileDir)
save(dirs.RecFileDir, 'DRF'); %Only save if it was a full set of trials

%Draw the OST progression, if you want to
if expParam.bVis == 1
    OST_MULT = 500; %Scale factor for OST
    visSignals(data, 16000, OST_MULT, savedWavdir)
end
end

function run = prompt4RunName()

prompt = 'Name of Run?:';
name   = 'Run Name';
numlines = 1;
defaultanswer = {'SF'};
runPrompt = inputdlg(prompt, name, numlines, defaultanswer);

run = runPrompt{1};
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

function dispLoudnessResult(loudResult)

switch loudResult
    case -1
        result = 'too soft';
    case 0
        result = 'just right';
    case 1
        result = 'too loud';
end

fprintf('Subject was %s\n', result)
end

function [f0b, targRMS, rmsB, gender, age] = loadBaselineVoice(dirs)

if exist(dirs.BaseFile, 'File')
    load(dirs.BaseFile, 'DRF')
    
    f0b     = DRF.qRes.meanf0;
    targRMS = DRF.qRes.meanRMS;
    rmsB    = DRF.expParam.rmsB;
    gender  = DRF.expParam.gender;
    age     = DRF.expParam.age;
else
    fprintf('Could not find baseline voice file at %s\n', dirs.BaseFile)
    fprintf('Loading Default Values for f0b, meanRMS, and rmsB\n')
    f0b     = 100;
    targRMS = 70.00;
    rmsB    = 0.00002;
    gender  = 'female';
    age     = 20;
end
end