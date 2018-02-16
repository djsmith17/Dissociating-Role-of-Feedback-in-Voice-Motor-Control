function dfRecBaselineVoice()
% dfRecBaselineVoice() is a simple script for recording short samples of a 
% speakers voice using Audapter and a MOTU Audio Card. The output of these 
% recordings sill give you an average RMS, and pitch of the samples. 
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

close all;

prompt = {'Subject ID:',...
          'Session ID:',...
          'Gender ("male" or "female")',...
          'Number of Trials:'};
name = 'Subject Information';
numlines = 1;
defaultanswer = {'null', 'BV1', 'female', '3'};
answer = inputdlg(prompt, name, numlines, defaultanswer);

if isempty(answer)
    return
end

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
expParam.subject    = answer{1}; %Subject#, Pilot#, null
expParam.run        = answer{2};
expParam.curSess    = [expParam.subject expParam.run];
expParam.gender     = answer{3};
expParam.numTrial   = str2double(answer{4});
expParam.AudFBSw    = 0;
expParam.trialLen   = 4; %Seconds

dirs = dfDirs(expParam.project);

dirs.RecFileDir = fullfile(dirs.RecData, expParam.subject, expParam.run);
dirs.RecWaveDir = fullfile(dirs.RecFileDir, 'wavFiles');

if exist(dirs.RecFileDir, 'dir') == 0
    mkdir(dirs.RecFileDir)
end
if exist(dirs.RecWaveDir, 'dir') == 0
    mkdir(dirs.RecWaveDir)
end

expParam.sRate              = 48000;  % Hardware sampling rate (before downsampling)
expParam.frameLen           = 96;      % Before downsampling
expParam.downFact           = 3;
expParam.sRateAnal          = expParam.sRate/expParam.downFact;
expParam.frameLenDown       = expParam.frameLen/expParam.downFact;
expParam.audioInterfaceName = 'MOTU MicroBook'; % 'ASIO4ALL' 'Komplete'

%Set up Audapter
Audapter('deviceName', expParam.audioInterfaceName);
Audapter('setParam', 'downFact', expParam.downFact, 0);
Audapter('setParam', 'sRate', expParam.sRateAnal, 0);
Audapter('setParam', 'frameLen', expParam.frameLen / expParam.downFact, 0);
p = getAudapterDefaultParams(expParam.gender);

%Set up OST and PCF Files.
expParam.ostFN = fullfile(dirs.Prelim, 'SFPerturbOST.ost'); check_file(expParam.ostFN);
expParam.pcfFN = fullfile(dirs.Prelim, 'SFPerturbPCF.pcf'); check_file(expParam.pcfFN);

[expParam, p]      = dfSetAudFB(expParam, dirs, p); %Sets some p params

expParam.boundsRMS = 3;
expParam.targRMS   = 60;

if VoiceRecsw == 1
    expParam.cuePause = 1.0;
    expParam.resPause = 2.0;
else
    expParam.cuePause = 0;
    expParam.resPause = 0;
end

%%%%%Visual Presentation
[~, H1, H2, H3, ~, ~, trigCirc] = dfSetVisFB(expParam.targRMS, expParam.boundsRMS);

%Open the curtains
pause(5); %Let them breathe a sec
set(H3,'Visible','off');

rawData = []; rmsData = [];
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
    set([H2 trigCirc],'Visible','on');
    
    fprintf('Trial %d\n',ii)
    AudapterIO('init', p);
    Audapter('reset');
    Audapter('start');

    pause(expParam.trialLen);
    
    Audapter('stop');
    set([H2 trigCirc],'Visible','off');
    
    % Save the data
    data    = dfSaveRawData(expParam, dirs);
    rawData = cat(1, rawData, data);
    
    pause(expParam.resPause)
end
close all

DRF.dirs        = dirs;
DRF.expParam    = expParam;
DRF.p           = p;
DRF.rawData     = rawData; 

%Do some quick analysis
qRes     = dfAnalysisAudioQuick(DRF, 1);
DRF.qRes = qRes;

dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.subject expParam.run dirs.saveFileSuffix 'DRF.mat']);
save(dirs.RecFileDir, 'DRF')

fprintf('\nThe mean f0 of each recordings were\n %4.2f Hz, %4.2f Hz, and %4.2f Hz\n', qRes.trialf0)
fprintf('\nThe mean f0 of all voice recordings\n is %4.2f Hz\n', qRes.meanf0)

fprintf('\nThe mean Amplitude of each recordings were\n %4.2f dB, %4.2f dB, and %4.2f dB\n', qRes.audioRMS)
fprintf('\nThe mean Amplitude of all voice recordings\n is %4.2f dB\n\n', qRes.meanRMS)
end