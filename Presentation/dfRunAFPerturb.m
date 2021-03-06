function dfRunAFPerturb()
% dfRunAFPerturb()
% Pitch-shift Perturbation experiment. This script measures acoustic output 
% from a participant as they have their auditory feedback perturbed.
% Audapter collects and %manages the recorded acoustic data. This 
% specifically uses a pitch-shift that matches the size of the stimulus seen
% in the somatosensory perturbation experiment. 
%
% This script calls the following 9 functions:
% dfDirs.m
% initNIDAQ.m
% dfSetAudFB.m
% dfSetTrialOrder.m
% dfMakePertSignal.m
% dfSetAudapFiles.m
% dfSetVisFB.m
% dfSaveRawData.m
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
lenDb = 0;
boxPos = setDialBoxPos(lenDb);

% Main Experimental prompt: Subject/Run Information
subject    = 'Pilot0';
run        = prompt4RunName();
InflaVarNm = 'IV1';
baseV      = 'BV1';

% Dialogue box asking for what type of Pitch-Shifted Feedback?
pertType = 'Linear Standard'; %questdlg('What type of Perturbation?', 'Type of Perturbation?', 'Linear Standard', 'Sinusoid Matched', 'Sinusoid Matched');
switch pertType
    case 'Linear Standard'
        pertTypeSw = 0;
    case 'Sinusoid Matched'
        pertTypeSw = 1;
end

AlgoType = 'pp_none'; %MFquestdlg(boxPos, 'What type of Perturbation?', 'Type of Perturbation?', 'pp_none', 'pp_peaks', 'pp_valleys', 'pp_none');

% Dialogue box asking if Practice set or Full set of trials
recType = MFquestdlg(boxPos, 'Practice or Full?','Length', 'Practice', 'Diagnostic', 'Full','Full');
switch recType
    case 'Practice'
        numTrials = 4;
        perCatch  = 1.00;
    case 'Diagnostic'
        numTrials = 10;
        perCatch  = 0.5;
    case 'Full'
        numTrials = 40;
        perCatch  = 0.25;
end

%Experiment Configurations
expParam.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType      = 'Auditory Perturbation_Perceptual';
expParam.subject      = subject;
expParam.run          = run;
expParam.curSess      = [expParam.subject expParam.run];
expParam.balloon      = 'N/A';
expParam.tightness    = 'N/A';
expParam.InflaVarNm   = InflaVarNm;
expParam.niDev        = 'Dev1';                      % NIDAQ Device Name. For more information, see dfInitNIDAQ
expParam.trialLen     = 4;                           % Seconds
expParam.numTrial     = numTrials;
expParam.curTrial     = [];
expParam.perCatch     = perCatch;
expParam.headGain     = 5;                   % Output gain above the input
expParam.AudFB        = 'Voice Shifted';
expParam.AudFBSw      = 1; %Voice Shifted
expParam.AudPert      = pertType;
expParam.AudPertSw    = pertTypeSw;
expParam.rmsThresh    = 0.011;
expParam.pitchShiftAlgo = AlgoType;
expParam.bVis         = 0;
expParam.sensorPType  = 'Seven';

%Set our dirs based on the project
dirs = dfDirs(expParam.project);
% Folder paths to save data files
dirs.RecFileDir = fullfile(dirs.RecData, expParam.subject, expParam.run);
dirs.RecWaveDir = fullfile(dirs.RecFileDir, 'wavFiles');
dirs.BaseFile   = fullfile(dirs.RecData, expParam.subject, baseV, [expParam.subject baseV 'DRF.mat']);

if exist(dirs.RecFileDir, 'dir') == 0
    mkdir(dirs.RecFileDir)
end
if exist(dirs.RecWaveDir, 'dir') == 0
    mkdir(dirs.RecWaveDir)
end

[expParam.f0b,...
 expParam.targRMS,...
 expParam.rmsB, ...
 expParam.gender,...
 expParam.age] = loadBaselineVoice(dirs);

% Look for the Inflation Response Files. Should Return InflaVar
expParam.InflaFile = [expParam.subject expParam.InflaVarNm 'DRF.mat']; % Results from the laryngeal perturbation experiment
dirs.InflaVarFile  = fullfile(dirs.SavData, expParam.subject, expParam.InflaVarNm, expParam.InflaFile);
if ~exist(dirs.InflaVarFile, 'file')
    fprintf('Warning: No Inflation Vars File at %s!\n', dirs.InflaVarFile)
    fprintf('Will use default Inflation Vars instead\n')
    InflaVar = [0.100 -100];
else
    fprintf('Inflation Variables found!!\n')
    load(dirs.InflaVarFile);
end

%Paradigm Configurations
expParam.sRate              = 48000;  % Hardware sampling rate (before downsampling)
expParam.frameLen           = 192;     % Before downsampling
expParam.downFact           = 3;
expParam.sRateAnal          = expParam.sRate/expParam.downFact; %Everything get automatically downsampled! So annoying
expParam.frameLenDown       = expParam.frameLen/expParam.downFact;
expParam.audioInterfaceName = 'MOTU MicroBook'; %'ASIO4ALL' 'Komplete'

%Set up Audapter
% Audapter('deviceName', expParam.audioInterfaceName);
% Audapter('setParam', 'downFact', expParam.downFact, 0);
% Audapter('setParam', 'sRate', expParam.sRateAnal, 0);
% Audapter('setParam', 'frameLen', expParam.frameLenDown, 0);
p = getAudapterDefaultParams(expParam.gender);
p.rmsThresh        = expParam.rmsThresh;
p.frameLen         = expParam.frameLenDown;
p.timeDomainPitchShiftAlgorithm = AlgoType;

%Set up Parameters to control NIDAQ and Perturbatron
[s, niCh, nVS]  = dfInitNIDAQ(expParam.niDev, expParam.trialLen);
expParam.sRateQ = s.Rate; % NIDAQ sampling rate
expParam.niCh   = niCh;   % Structure of Channel Names

% %Set up OST and PCF Files
% expParam.ostFN = fullfile(dirs.Prelim, 'AFPerturbOST.ost'); check_file(expParam.ostFN);
% expParam.pcfFN = fullfile(dirs.Prelim, 'AFPerturbPCF.pcf'); check_file(expParam.pcfFN);

%Set up Auditory Feedback (Masking Noise, Pitch-Shift?)
[p, SSNw, SSNfs]  = dfSetAudFB(expParam, dirs, p);

%Set up the order of trials (Order of perturbed, control, etc)
expParam.trialType = dfSetTrialOrder(expParam.numTrial, expParam.perCatch);

%Select the trigger points for perturbation onset and offset and creating
%the digital signal to be sent to the NIDAQ
[expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, lenDb);

expParam.rdyPause  = 5.0; % How long to show them 'Ready'
expParam.cuePause  = 1.0; % How long the cue period lasts
expParam.buffPause = 0.8; % Give them a moment to start speaking
expParam.endPause  = 0.5;
expParam.resPause  = 2.0; % How long the rest/VisFB lasts
expParam.boundsRMS = 3;   % +/- dB

% Load the InflaVar Variables. Analyzed from previous recording
expParam.InflaT   = InflaVar(1);
expParam.InflaV   = InflaVar(2);

%This is where the fun begins
fprintf('\nStarting Trials\n\n')

%Dim the lights (Set the visual Feedback)
[msrStr, annoStr] = dfSetVisFB(2, expParam.curSess, expParam.targRMS, expParam.boundsRMS);

% Only play masking noise for this condition
if expParam.AudFBSw == 2
    sound(SSNw, SSNfs)
end

%Open the curtains
pause(expParam.rdyPause);            % Let them breathe a sec
set(annoStr.Ready, 'Visible','off'); % Turn off 'Ready?'

DAQin = []; rawData = [];
loudResults = [];
for ii = 1:expParam.numTrial
    expParam.curTrialNum  = ii;
    expParam.curTrial     = ['Trial' num2str(ii)];
    expParam.curSessTrial = [expParam.subject expParam.run expParam.curTrial];
    
    %Level of f0 change based on results from Laryngeal pert Exp
    audStimP = dfSetAudapFiles(dirs, expParam, ii);
    p.timeDomainPitchShiftSchedule = audStimP.pertSched; 
    
    %Set the OST and PCF functions
%     Audapter('ost', expParam.ostFN, 0);
%     Audapter('pcf', expParam.pcfFN, 0);
    
    %Setup which perturb file we want
    NIDAQsig = [expParam.sigs(:,ii) nVS];
    queueOutputData(s, NIDAQsig);
    
    %Cue to begin trial
    set(annoStr.plus, 'Visible','on');
    pause(expParam.cuePause)
    
    %Phonation Start
    set(annoStr.plus, 'Visible','off');
    set([annoStr.EEE annoStr.visTrig],'Visible','on');
   
    fprintf('Trial %d\n',ii)
    AudapterIO('init', p);
    Audapter('reset');
    Audapter('start');
    Audapter('stop');
    Audapter('start');
    
%     pause(expParam.buffPause)

    %Play out the Analog Perturbatron Signal. This will hold script for as
    %long as vector lasts. In this case, 4.0 seconds. 
    [dataDAQ, ~] = s.startForeground;
    
    %Phonation End
    set([annoStr.EEE annoStr.visTrig],'Visible','off');
    pause(expParam.endPause)
    Audapter('stop');
    
    %Save the data
    data    = AudapterIO('getData'); % This will need to become a try statement again
    DAQin   = cat(3, DAQin, dataDAQ);
    rawData = cat(1, rawData, data);
    
    %Grab smooth RMS trace from 'data' structure
    rmsMean = dfCalcMeanRMS(data, expParam.rmsB);
    %Compare against baseline and updated Visual Feedback
    [color, newPos, loudResult] = dfUpdateVisFB(msrStr, rmsMean);
    loudResults = cat(1, loudResults, loudResult);
    dispLoudnessResult(loudResult)
    
    set(annoStr.LoudRec, 'position', newPos);
    set(annoStr.LoudRec, 'Color', color, 'FaceColor', color);
    set([annoStr.LoudRec annoStr.fbLines], 'Visible', 'on'); 

    dfSaveWavRec(data, expParam, dirs);
    pause(expParam.resPause)
    set([annoStr.LoudRec annoStr.fbLines], 'Visible', 'off');
end
close all;
elapsed_time = toc(ET)/60;    % Elapsed Time of the experiment
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
        DRF.qRes = dfAnalysisAudioQuick(DRF, 1);
end

% Save the large structure (only if not practice trials)
dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.subject expParam.run dirs.saveFileSuffix 'DRF.mat']);
fprintf('\nSaving recorded data at:\n%s\n\n', dirs.RecFileDir)
save(dirs.RecFileDir, 'DRF'); %Only save if it was a full set of trials

dirs.SavResultsDir = fullfile(dirs.Results, expParam.subject, expParam.run);
if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

f0b = 100;
aFa = 1; iRf = 0;
niAn = struct;
niAn.sRate = 8000;
% [~, auRes] = dfAnalysisAudapter(dirs, DRF.expParam, DRF.rawData, f0b, aFa, iRf, niAn);
% 
% drawAudRespMeanTrial(auRes, dirs.SavResultsDir)
% pause(2)
% drawAudRespIndivTrial(auRes, dirs.SavResultsDir)
% pause(2)

%Draw the OST progression, if you want to
if expParam.bVis == 1
    OST_MULT = 500; %Scale factor for OST
    visSignals(data, 16000, OST_MULT, savedWavdir)
end
end

function boxPos = setDialBoxPos(debug)

monitorSize = get(0, 'Monitor');
numMon = size(monitorSize, 1);

boxPos = [0.45 0.45];

if debug == 1 && numMon > 1
    boxPos = boxPos + [1 0];
end
end

function run = prompt4RunName()

prompt = 'Name of Run?:';
name   = 'Run Name';
numlines = 1;
defaultanswer = {'AF'};
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