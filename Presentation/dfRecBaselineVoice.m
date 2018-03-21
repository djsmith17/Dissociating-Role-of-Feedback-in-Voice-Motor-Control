function dfRecBaselineVoice()
% dfRecBaselineVoice() is a simple script for recording short samples of a 
% speakers voice using Audapter and a MOTU Audio Card. The output of these 
% recordings will give you an average RMS, and pitch of the samples. 
% This should be used at the beginning of an audapter recording session to 
% determine a baseline voice amplitude and pitch. 
% 
% This can also be used to calibrate the microphone. There is a switch case
% that can be activated when the script is run. Make sure to save it as MC
% instead of BV.
%
% This script assumes:
% 1: The participant speaks at a comfortable and typical speaking volume
% 2: The microphone is placed at a fixed distance (e.g. 7cm) from the participant
% 3: The microphone gain levels are constant for each participant and through the trials
% 4: The participant phonates a steady-state vowel sound through these recordings
%
% This script calls the following 5 functions:
% dfDirs.m
% dfSetAudFB.m
% dfSetVisFB.m
% dfSaveRawData.m
% dfAnalysisAudioQuick.m
%
% This uses the toolbox from MATLAB-Toolboxes
% speechres
%
% This script is also dependent on the following Mathworks Toolboxes
% Signal-Processing Toolbox

close all;

% Main Experimental prompt: Subject/Run Information
subject    = 'Pilot33'; % Subject#, Pilot#, null
run        = 'BV1';     % Baseline Voice (BV) or Calibrate Microphone (CM)
gender     = 'male';    % "male" or "female"
numTrials  = 3;         % number of trials;

VoiceRec = questdlg('Calibrate Mic or Baseline Voice?', 'Recording Type', 'Calibrate Microphone', 'Baseline Voice', 'Baseline Voice');
switch VoiceRec
    case 'Calibrate Microphone'
        VoiceRecsw = 0;
    case 'Baseline Voice'
        VoiceRecsw = 1;
end

%Paradigm Configurations
expParam.project    = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType    = 'Somatosensory Perturbation_Perceptual';
expParam.subject    = subject; 
expParam.run        = run;
expParam.curSess    = [expParam.subject expParam.run];
expParam.gender     = gender;

if VoiceRecsw == 1
    expParam.trialLen = 4;                        % Seconds
    expParam.numTrial = numTrials;
    expParam.AudFBSw  = 0;
    expParam.cuePause = 1.0;
    expParam.resPause = 2.0;
else
    expParam.trialLen = 100;                      % Seconds
    expParam.numTrial = 1;
    expParam.AudFBSw  = 1;
    expParam.cuePause = 0;
    expParam.resPause = 0;
end

fprintf('\nBeginning baseline voice recordings for\n%s %s\n\n', subject, run)

% Set our dirs based on the project
dirs = dfDirs(expParam.project);

dirs.RecFileDir = fullfile(dirs.RecData, expParam.subject, expParam.run);
dirs.RecWaveDir = fullfile(dirs.RecFileDir, 'wavFiles');

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
expParam.sRateAnal          = expParam.sRate/expParam.downFact;
expParam.frameLenDown       = expParam.frameLen/expParam.downFact;
expParam.audioInterfaceName = 'MOTU MicroBook'; % 'ASIO4ALL' 'Komplete'

%Set up Audapter
Audapter('deviceName', expParam.audioInterfaceName);
Audapter('setParam', 'downFact', expParam.downFact, 0);
Audapter('setParam', 'sRate', expParam.sRateAnal, 0);
Audapter('setParam', 'frameLen', expParam.frameLenDown, 0);
p = getAudapterDefaultParams(expParam.gender);

%Set up OST and PCF Files.
expParam.ostFN = fullfile(dirs.Prelim, 'SFPerturbOST.ost'); check_file(expParam.ostFN);
expParam.pcfFN = fullfile(dirs.Prelim, 'SFPerturbPCF.pcf'); check_file(expParam.pcfFN);

% Set up Auditory Feedback (Voice Not Shifted)
[expParam, p]      = dfSetAudFB(expParam, dirs, p);

expParam.boundsRMS = 3;
expParam.targRMS   = 70;

% Dim the lights (Set the visual Feedback)
[~, H1, H2, H3, ~, ~, trigCirc] = dfSetVisFB(expParam.curSess, expParam.targRMS, expParam.boundsRMS);

%Open the curtains
pause(5);                % Let them breathe a sec
set(H3,'Visible','off'); % Turn off 'Ready?'

rawData = [];
for ii = 1:expParam.numTrial
    expParam.curTrial     = ['Trial' num2str(ii)];
    expParam.curSessTrial = [expParam.subject expParam.run expParam.curTrial];
    
    %Set the OST and PCF functions
    Audapter('ost', expParam.ostFN, 0);
    Audapter('pcf', expParam.pcfFN, 0);
    
    %Cue to begin trial
    set(H1,'Visible','on');
    pause(expParam.cuePause)
    
    %Phonation Start
    set(H1,'Visible','off');
    set([H2 trigCirc],'Visible','on'); % Turn on the 'eee' and trigMark
    
    fprintf('Trial %d\n',ii)
    AudapterIO('init', p);
    Audapter('reset');
    Audapter('start');

    pause(expParam.trialLen);
    
    Audapter('stop');
    set([H2 trigCirc],'Visible','off'); % Turn off the 'eee' and trigMark
    
    % Load the Audapter saved data and save some as wav Files
    data    = dfSaveRawData(expParam, dirs);
    rawData = cat(1, rawData, data);
    
    pause(expParam.resPause)
end
close all

% Store all the variables and data from the session in a large structure
DRF.dirs        = dirs;
DRF.expParam    = expParam;
DRF.p           = p;
DRF.rawData     = rawData; 

% Do some quick analysis
qRes     = dfAnalysisAudioQuick(DRF, 1);
DRF.qRes = qRes;

% Save the large data structure
dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.subject expParam.run dirs.saveFileSuffix 'DRF.mat']);
switch VoiceRec
    case 'Calibrate Microphone'
        return
    case 'Baseline Voice'
        fprintf('\nSaving recorded baseline data at:\n%s\n\n', dirs.RecFileDir)
        save(dirs.RecFileDir, 'DRF')
end

fprintf('\nThe mean f0 of each recordings were\n %4.2f Hz, %4.2f Hz, and %4.2f Hz\n', qRes.trialf0)
fprintf('\nThe mean f0 of all voice recordings\n is %4.2f Hz\n', qRes.meanf0)

fprintf('\nThe mean Amplitude of each recordings were\n %4.2f dB, %4.2f dB, and %4.2f dB\n', qRes.audioRMS)
fprintf('\nThe mean Amplitude of all voice recordings\n is %4.2f dB\n\n', qRes.meanRMS)
end