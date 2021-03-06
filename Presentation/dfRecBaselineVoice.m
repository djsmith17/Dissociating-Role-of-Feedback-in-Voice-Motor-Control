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
subject    = 'DRF12'; % Subject#, Pilot#, null
run        = 'BV1';
gender     = 'female';    % "male" or "female"
DOB        = datetime(1997, 1, 5); % Year, Month, Day
numTrials  = 3;         % number of trials;

% Paradigm Configurations
expParam.project    = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType    = 'Somatosensory Perturbation_Perceptual';
expParam.curDT      = datetime('now'); % Current Date and Time

% Set our dirs based on the project
dirs = dfDirs(expParam.project);

boxPos  = dfSetDialBoxPos(1);
recType = MFquestdlg(boxPos, 'Calibrate Mic or Baseline Voice?', 'Recording Type', 'Calibrate Microphone', 'Baseline Voice', 'Baseline Voice');
switch recType
    case 'Baseline Voice'
        expParam.subject  = subject; 
        expParam.run      = run;
        expParam.curSess  = [expParam.subject expParam.run];
        expParam.gender   = gender;
        expParam.DOB      = DOB;
        expParam.age      = years(expParam.curDT - expParam.DOB);
        expParam.f0b      = 150;
        expParam.trialLen = 4;                      % Seconds
        expParam.numTrial = numTrials;
        expParam.AudFBSw  = 0;
        expParam.rdyPause = 5;
        expParam.cuePause = 1.0;
        expParam.resPause = 2.0;
        expParam.headGain = 5;
        defMon = 2;
        
        expParam.rmsB     = loadCalibration(dirs, expParam.curDT);

        fprintf('\nBeginning baseline voice recordings for\n%s %s\n\n', subject, run)
    case 'Calibrate Microphone'
        expParam.subject  = 'Microphone Calibration';
        expParam.run      = 'MC';
        expParam.curSess  = expParam.run; 
        expParam.gender   = gender;
        expParam.DOB      = DOB;
        expParam.age      = years(expParam.curDT - expParam.DOB);
        expParam.f0b      = 150;
        expParam.trialLen = 34;                     % Seconds
        expParam.numTrial = 1;
        expParam.AudFBSw  = 0;
        expParam.rdyPause = 2;
        expParam.cuePause = 0;
        expParam.resPause = 0;
        expParam.headGain = 5;        
        defMon = 1;

        expParam.targLoud   = prompt4Calibrate(); % Ask user how loud the calibration sound source was
        expParam.measAmpInc = [];
    case ''
        fprintf('\nExperiment Exited\n')
        return
end

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
[p, ~, ~] = dfSetAudFB(expParam, dirs, p);

expParam.boundsRMS = 3;
expParam.targRMS   = 70;

% Dim the lights (Set the visual Feedback)
[~, annoStr] = dfSetVisFB(defMon, expParam.curSess, expParam.targRMS, expParam.boundsRMS);

%Open the curtains
pause(expParam.rdyPause); % Let them breathe a sec
set(annoStr.Ready,'Visible','off');  % Turn off 'Ready?'

rawData = [];
contTrials = 1; ii = 1;
while contTrials ~= 0
    
    expParam.curTrial     = ['Trial' num2str(ii)];
    expParam.curSessTrial = [expParam.curSess expParam.curTrial];
    
    %Set the OST and PCF functions
    Audapter('ost', expParam.ostFN, 0);
    Audapter('pcf', expParam.pcfFN, 0);
    
    %Cue to begin trial
    set(annoStr.plus, 'Visible','on');
    pause(expParam.cuePause)
    
    %Phonation Start
    set(annoStr.plus, 'Visible','off');
    set([annoStr.EEE annoStr.visTrig],'Visible','on'); % Turn on the 'eee' and visTrig
    
    fprintf('Trial %d\n',ii)
    AudapterIO('init', p);
    Audapter('reset');
    Audapter('start');

    pause(expParam.trialLen);
    
    Audapter('stop');
    set([annoStr.EEE annoStr.visTrig],'Visible','off'); % Turn off the 'eee' and visTrig
    
    % Load the Audapter saved data and save some as wav Files
    data = AudapterIO('getData');       % This will need to become a try statement again
    rawData = cat(1, rawData, data);
    
    dfSaveWavRec(data, expParam, dirs);
    
    pause(expParam.resPause)
    
    switch recType
        case 'Calibrate Microphone'
            rerun = MFquestdlg(boxPos, 'Rerun calibration step?', 'Rerun?', 'Yes', 'No', 'No');
            switch rerun
                case 'Yes'
                    rawData = [];
                case 'No'
                    contTrials = 0;
                    
                    expParam.measAmpInc = prompt4AmplifiedHead();
                    fprintf('The input microphone signal was measured at %0.2f dB\n', expParam.targLoud)
                    fprintf('The output headphone signal was measured at %0.2f dB\n', expParam.measAmpInc)
            end
        case 'Baseline Voice'
            ii = ii + 1;
            if ii > expParam.numTrial
                contTrials = 0;
            end
    end  
end
close all

% Store all the variables and data from the session in a large structure
DRF.dirs        = dirs;
DRF.expParam    = expParam;
DRF.p           = p;
DRF.rawData     = rawData; 

switch recType
    case 'Calibrate Microphone'
        DRF.expParam.rmsB = findRMS2dBRatio(rawData, expParam.targLoud);
    case 'Baseline Voice'
        DRF.qRes = dfAnalysisAudioQuick(DRF, 1); % Do some quick analysis
        displayBaselineResults(DRF.qRes)
end

% Save the large data structure
dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.curSess dirs.saveFileSuffix 'DRF.mat']);
fprintf('\nSaving recorded baseline data at:\n%s\n\n', dirs.RecFileDir)
save(dirs.RecFileDir, 'DRF')

dfInspectRawData(dirs.RecFileDir)
end

function targLoud = prompt4Calibrate()

prompt = 'Source Loudness (dB (SPL HL)):';
name = 'Source Loudness';
numlines = 1;
defaultanswer = {'75.00'};
loudnessPrompt = inputdlg(prompt, name, numlines, defaultanswer);

targLoud = str2double(loudnessPrompt{1});
end

function measAmp = prompt4AmplifiedHead()

prompt = 'Measured Insert Earphone Loudness (dB (SPL HL)):';
name = 'Earphone Loudness';
numlines = 1;
defaultanswer = {'80.00'};
loudnessPrompt = inputdlg(prompt, name, numlines, defaultanswer);

measAmp = str2double(loudnessPrompt{1});
end

function rmsB = loadCalibration(dirs, curDT)

%It will always be called MC/MCDRF.mat. However, we will check the age of
%the calibration
MicCalibFile = fullfile(dirs.RecData, 'Microphone Calibration', 'MC', 'MCDRF.mat');
if exist(MicCalibFile, 'file')
    fprintf('Loading the Microphone Calibration File.\n')
    load(MicCalibFile, 'DRF')

    rmsB  = DRF.expParam.rmsB;
    oldDT = DRF.expParam.curDT;
    dateCH = hours(curDT - oldDT);
else
    fprintf('The Microphone calibration file could not be found.\n')
    fprintf('Using default values for rmsB.\n')
    rmsB   = 0.00000002;
    dateCH = 0;
end

fprintf('Microphone Calibration is %0.1f hours old,\n', dateCH)
if dateCH > 24
    fprintf('You may want to retake the calibration\n')
end
fprintf('\n')
end

function displayBaselineResults(qRes)
fprintf('\nThe mean f0 of each recordings were\n %4.2f Hz, %4.2f Hz, and %4.2f Hz\n', qRes.trialf0)
fprintf('\nThe mean f0 of all voice recordings\n is %4.2f Hz\n', qRes.meanf0)

fprintf('\nThe mean Amplitude of each recordings were\n %4.2f dB, %4.2f dB, and %4.2f dB\n', qRes.audioRMS)
fprintf('\nThe mean Amplitude of all voice recordings\n is %4.2f dB\n\n', qRes.meanRMS)
end