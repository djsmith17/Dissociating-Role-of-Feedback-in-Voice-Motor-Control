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
subject    = 'Test';   % Subject#, Pilot#, null
run        = 'SFL1';
jitt       = 1;        % 0 Starts with Masking Noise, 1 Starts with Voice

balloon    = 'GB2';    % Which perturbation balloon?
baseV      = 'BVEndo';
FBNames    = {'Voice Feedback'; 'AC Masking Noise'};
FBTypes    = [0 2];
FBInstr    = {'Your Own Voice'; 'Loud White Noise'};

expParam = dfInitExpParam();

expParam.expType = 'Somatosensory Perturbation_Endoscopy';
expParam.subject   = subject;
expParam.run       = run;
expParam.curSess   = [expParam.subject expParam.run];
expParam.balloon   = balloon;
expParam.niDev     = 'Dev1';
expParam.numTrial  = 10; 
expParam.perCatch   = 1; % 100% of Trials
expParam.numMaskRep = 1;

expParam.rdyPause   = 2;

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
[msrStr, annoStr] = dfSetVisFB(2, expParam.curSess, expParam.targRMS, expParam.boundsRMS);

% Open the curtains
set([annoStr.curSessNote, annoStr.Ready], 'Visible', 'off') % Turn off 'Ready?' 

LR = LiveSensorResult(expParam, 1);
for ii = 1:expParam.numTrial
    ET = tic;
    
    % Do some Setup
    expParam.curTrialNum  = ii;
    expParam.curTrial     = ['Trial' num2str(ii)];
    expParam.curSessTrial = [expParam.subject expParam.run expParam.curTrial];    
    
    curStimType      = 2;
    expParam.AudFB   = FBNames{curStimType}; % Alternating trials
    expParam.AudFBSw = FBTypes(curStimType); % Alternating trials
    instrFB          = FBInstr(curStimType); % Alternating trials
    
    % Set up Auditory Feedback (Masking Noise, Pitch-Shift?)
    [p, SSNw, SSNfs]      = dfSetAudFB(expParam, dirs, p);    
        
    %Set the OST and PCF functions
    Audapter('ost', expParam.ostFN, 0);
    Audapter('pcf', expParam.pcfFN, 0);
    
    %Setup which perturb file we want
    NIDAQsig = [expParam.sigs(:,ii) nVS];
    queueOutputData(s, NIDAQsig);
    
    %Pause and get ready to give them the next trial
    DAQin       = []; rawData  = [];
    loudResults = []; audStimP = [];
    
    set(annoStr.trialT, 'String', ['Trial ' num2str(ii) '/' num2str(expParam.numTrial)]);
    set(annoStr.FBCue, 'String', instrFB);
    set([annoStr.trialT annoStr.FBCue], 'Visible', 'on');
    pause()
    set([annoStr.trialT annoStr.FBCue], 'Visible', 'off');
    
    % Only play masking noise for this condition
    if expParam.AudFBSw == 2
        sound(SSNw, SSNfs)
    end
    pause(expParam.rdyPause)
    
    %Cue to begin trial
    set(annoStr.plus,'Visible','on');
    pause(expParam.cuePause)
    
    %Phonation Start
    set(annoStr.plus,'Visible','off');
    set(annoStr.EEE,'Visible','on');
    
    fprintf('Trial %d\n',ii)
    AudapterIO('init', p);
    Audapter('reset');
    Audapter('start');
    pause(expParam.buffPause)
    
    % Play out the Analog Perturbatron Signal. This will hold script for as
    % long as vector lasts. In this case, 4.0 seconds. 
    [dataDAQ, ~] = s.startForeground;
     
    % Phonation End
    set(annoStr.EEE,'Visible','off');
    pause(expParam.endPause)
    Audapter('stop');
    
    % Load the Audapter saved data and save as wav Files
    data    = AudapterIO('getData'); % This will need to become a try statement again
    DAQin   = cat(3, DAQin, dataDAQ);
    rawData = cat(1, rawData, data);    
       
    % Grab smooth RMS trace from 'data' structure
    rmsMean = dfCalcMeanRMS(data, expParam.rmsB);
    % Compare against baseline and updated Visual Feedback
    [color, newPos, loudResult] = dfUpdateVisFB(msrStr, rmsMean);
    loudResults = cat(1, loudResults, loudResult);
    dispLoudnessResult(loudResult)

    set(annoStr.LoudRec, 'position', newPos);
    set(annoStr.LoudRec, 'Color', color, 'FaceColor', color);
    set([annoStr.LoudRec annoStr.fbLines], 'Visible', 'on');   
    pause(expParam.resPause)
    set([annoStr.LoudRec annoStr.fbLines], 'Visible', 'off');
    
    LR = LR.updateLiveResult(dataDAQ, ii);    
    dfSaveWavRec(data, expParam, dirs); 
    
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