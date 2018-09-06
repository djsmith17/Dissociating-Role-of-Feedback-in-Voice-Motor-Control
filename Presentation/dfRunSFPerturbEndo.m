function dfRunSFPerturbEndo()
% dfRunSFPerturbEndo()
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
rng('shuffle');

% Main Experimental prompt: Subject/Run Information
subject    = 'DRF_EN0';   % Subject#, Pilot#, null
run        = 'SFL1';

balloon    = '2E4';     % Which perturbation balloon?
baseV      = 'BVEndo';
FBNames = {'Voice Feedback'; 'AC Masking Noise'};
FBTypes = [0 2];

expParam = dfInitExpParam();

expParam.subject   = subject;
expParam.run       = run;
expParam.curSess   = [expParam.subject expParam.run];
expParam.balloon   = balloon;
expParam.numTrial  = 2; 
expParam.perCatch  = 1; % 100% of Trials

% Set our dirs based on the project
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

[expParam.gender,...
 expParam.age,...
 expParam.f0b,...
 expParam.targRMS,...
 expParam.rmsB]         = loadBaselineVoice(dirs);

% Set up Audapter Parameters and function paths
[expParam, p] = dfInitAudapater(dirs, expParam);

% Set up NIDAQ Parameters and channel specifics
[s, niCh, nVS]  = dfInitNIDAQ(expParam.niDev, expParam.trialLen);
expParam.sRateQ = s.Rate; % NIDAQ sampling rate
expParam.niCh   = niCh;   % Structure of Channel Names

% Set up the order of trials (Order of perturbed, control, etc)
expParam.trialType = dfSetTrialOrder(expParam.numTrial, expParam.perCatch);

% Select the trigger points for perturbation onset and offset and creating
% the digital signal to be sent to the NIDAQ
[expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType);

% This is where the fun begins
fprintf('\nStarting Trials\n\n')

% Dim the lights (Set the visual Feedback)
[anMsr, H1, H2, H3, fbLines, rec, trigCirc] = dfSetVisFB(expParam.curSess, expParam.targRMS, expParam.boundsRMS);

% Open the curtains
pause(expParam.rdyPause);  % Let them breathe a sec
set(H3, 'Visible', 'off'); % Turn off 'Ready?'

LR = LiveSensorResult(expParam, 1);
for ii = 1:expParam.numTrial
    ET = tic;
    
    % Do some Setup
    expParam.curTrialNum  = ii;
    expParam.curTrial     = ['Trial' num2str(ii)];
    expParam.curSessTrial = [expParam.subject expParam.run expParam.curTrial];    
    
    expParam.AudFB   = FBNames(mod(ii,2)+ 1); % Alternating trials
    expParam.AudFBSw = FBTypes(mod(ii,2)+ 1); % Alternating trials
    
    % Set up Auditory Feedback (Masking Noise, Pitch-Shift?)
    [expParam, p]      = dfSetAudFB(expParam, dirs, p);    
        
    %Set the OST and PCF functions
    Audapter('ost', expParam.ostFN, 0);
    Audapter('pcf', expParam.pcfFN, 0);
    
    %Setup which perturb file we want
    NIDAQsig = [expParam.sigs(:,ii) nVS];
    queueOutputData(s, NIDAQsig);
    
    %Pause and get ready to give them the next trial
    DAQin       = []; rawData  = [];
    loudResults = []; audStimP = [];
    pause()
    
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
    
    % Play out the Analog Perturbatron Signal. This will hold script for as
    % long as vector lasts. In this case, 4.0 seconds. 
    [dataDAQ, ~] = s.startForeground;
     
    % Phonation End
    set([H2 trigCirc],'Visible','off');
    pause(expParam.endPause)
    Audapter('stop');
    
    % Load the Audapter saved data and save as wav Files
    data    = AudapterIO('getData'); % This will need to become a try statement again
    DAQin   = cat(3, DAQin, dataDAQ);
    rawData = cat(1, rawData, data);    
       
    % Grab smooth RMS trace from 'data' structure
    rmsMean = dfCalcMeanRMS(data, expParam.rmsB);
    % Compare against baseline and updated Visual Feedback
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
    
    elapsed_time = toc(ET)/60;   % Elapsed Time of the trial
    
    % Store all the variables and data from the session in a large structure
    expParam.elapsedTime = elapsed_time;
    expParam.loudResults = loudResults;
    DRF.dirs        = dirs;
    DRF.expParam    = expParam;
    DRF.p           = p;
    DRF.audStimP    = audStimP;
    DRF.DAQin       = DAQin;
    DRF.rawData     = rawData;     
    
    % Save the large data structure (only if not practice trials)
    dirs.RecFile = fullfile(dirs.RecFileDir, [expParam.subject expParam.run expParam.curTrial dirs.saveFileSuffix 'DRF.mat']);
    fprintf('\nSaving recorded data at:\n%s\n\n', dirs.RecFile)
    save(dirs.RecFile, 'DRF');    
end
close all;
end

function run = prompt4RunName()

prompt = 'Name of Run?:';
name   = 'Run Name';
numlines = 1;
defaultanswer = {'SF'};
runPrompt = inputdlg(prompt, name, numlines, defaultanswer);

run = runPrompt{1};
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

function [gender, age, f0b, targRMS, rmsB] = loadBaselineVoice(dirs)

if exist(dirs.BaseFile, 'File')
    load(dirs.BaseFile, 'DRF')
    
    gender  = DRF.expParam.gender;
    age     = DRF.expParam.age;
    f0b     = DRF.qRes.meanf0;
    targRMS = DRF.qRes.meanRMS;
    rmsB    = DRF.expParam.rmsB;
else
    fprintf('Could not find baseline voice file at %s\n', dirs.BaseFile)
    fprintf('Loading Default Values for f0b, meanRMS, and rmsB\n')
    
    gender  = 'female';
    age     = 20;
    f0b     = 200;
    targRMS = 70.00;
    rmsB    = 0.00002;
end
end