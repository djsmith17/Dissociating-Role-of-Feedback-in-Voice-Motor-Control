function dfRunAFPerturbOfflineTD()
% dfRunAFPerturbOffline()
% This script pitch-shifts a previously recorded wav file, using the same
% Audapter methods used during online trials. This script acts as a test 
% bed to make sure that everything works the way it should. 
%
% This script calls the following (8) functions:
% dfDirs.m
% dfInitNIDAQ.m
% dfSetAudFB.m
% dfSetTrialOrder.m
% dfMakePertSignal.m
% dfSetAudapFiles.m
% dfSaveRawData
%
% This uses the toolbox from MATLAB-Toolboxes
% speechres
%
% This script is also dependent on the following Mathworks Toolboxes
% Signal-Processing Toolbox

close all;
ET = tic;
rng('shuffle');
lenDb = 1;
boxPos = setDialBoxPos(lenDb);

% Main Experimental prompt: Subject/Run Information
subject    = 'Pilot0';
run        = prompt4RunName();
InflaVarNm = 'IV1';
BaseRun    = 'BV1';

LoadSavDataLoc = 0;
collectNewData = 1; %Boolean

% Dialogue box asking for what type of Pitch-Shifted Feedback?
pertType = 'Linear Standard'; %questdlg('What type of Perturbation?', 'Type of Perturbation?', 'Linear Standard', 'Sigmoid Matched', 'Sigmoid Matched');
switch pertType
    case 'Linear Standard'
        pertTypeSw = 0;
    case 'Sigmoid Matched'
        pertTypeSw = 1;
end

AlgoType = MFquestdlg(boxPos, 'What type of Perturbation?', 'Type of Perturbation?', 'pp_none', 'pp_peaks', 'pp_valleys', 'pp_none');

% Dialogue box asking if Practice set or Full set of trials
recType = MFquestdlg(boxPos, 'Practice or Full?','Length', 'Practice', 'Diagnostic', 'Full','Full');
switch recType
    case 'Practice'
        numTrials = 4;
        perCatch  = 1.0;
    case 'Diagnostic'
        numTrials = 10;
        perCatch  = 0.50;
    case 'Full'
        numTrials = 40;
        perCatch  = 0.25;
end

%Experiment Configurations
expParam.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType      = 'Auditory Perturbation_Perceptual';
expParam.subject      = subject;
expParam.run          = [run 'Offline'];
expParam.curSess      = [expParam.subject expParam.run];
expParam.balloon      = 'N/A';
expParam.tightness    = 'N/A';
expParam.InflaVarNm   = InflaVarNm;
expParam.niDev        = 'Dev2';          % NIDAQ Device Name. For more information, see dfInitNIDAQ
expParam.trialLen     = 4;               % Seconds
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
expParam.bVis         = 1;
expParam.bPlay        = 0;

%Set our dirs based on the project
dirs = dfDirs(expParam.project);
% Folder paths to save data files
dirs.RecFileDir = fullfile(dirs.RecData, expParam.subject, expParam.run);
dirs.RecWaveDir = fullfile(dirs.RecFileDir, 'wavFiles');
dirs.BaseFile   = fullfile(dirs.RecData, expParam.subject, BaseRun, [expParam.subject BaseRun 'DRF.mat']);

if exist(dirs.RecFileDir, 'dir') == 0
    mkdir(dirs.RecFileDir)
end
if exist(dirs.RecWaveDir, 'dir') == 0
    mkdir(dirs.RecWaveDir)
end

[expParam.f0b,...
 expParam.targRMS,...
 expParam.rmsB, ...
 expParam.gender] = loadBaselineVoice(dirs);

% Look for the Inflation Response Files
expParam.InflaFile = [expParam.subject expParam.InflaVarNm 'DRF.mat'];
dirs.InflaVarFile  = fullfile(dirs.RecData, expParam.subject, expParam.InflaVarNm, expParam.InflaFile);
if ~exist(dirs.InflaVarFile, 'file')
    fprintf('Warning: No Inflation Vars File at %s!\n', dirs.InflaVarFile)
    fprintf('Will use default Inflation Vars instead\n')
    InflaVar = [0.100 -100];
else
    fprintf('Inflation Variables found!!\n')
    load(dirs.InflaVarFile);
end

% Look for a place to save the data
dirs.SavResultsDir = fullfile(dirs.Results, expParam.subject, expParam.run);
if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

if collectNewData == 1
    %Paradigm Configurations
    expParam.sRate              = 48000;  % Hardware sampling rate (before downsampling)
    expParam.frameLen           = 192;  % Before downsampling
    expParam.downFact           = 3;
    expParam.sRateAnal          = expParam.sRate/expParam.downFact; %Everything get automatically downsampled! So annoying
    expParam.frameLenDown       = expParam.frameLen/expParam.downFact;
    expParam.audioInterfaceName = 'MOTU MicroBook'; %'ASIO4ALL' 'Komplete'

    %Set up Parameters to control NIDAQ and Perturbatron
    expParam.sRateQ = 8000; % NIDAQ sampling rate
    expParam.niCh   = [];   % Structure of Channel Names
    
    %Set up Audapter
    p = getAudapterDefaultParams(expParam.gender);
    p.rmsThresh        = expParam.rmsThresh;
    p.frameLen         = expParam.frameLenDown;
    
    %Set up Auditory Feedback (Masking Noise, Pitch-Shift?)
    [expParam, p]      = dfSetAudFB(expParam, dirs, p);    
    
    %Set up the order of trials (Order of perturbed, control, etc)
    expParam.trialType = dfSetTrialOrder(expParam.numTrial, expParam.perCatch); %numTrials, percentCatch

    %Select the trigger points for perturbation onset and offset and creating
    %the digital signal to be sent to the NIDAQ
    [expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, lenDb);

    expParam.cuePause  = 1.0; % How long the cue period lasts
    expParam.buffPause = 0.2; % Give them a moment to start speaking
    expParam.resPause  = 2.0; % How long the rest/VisFB lasts
    expParam.boundsRMS = 3;   % +/- dB
    
    % Gives variable of InflaVar. Analyzed from previous recording
    expParam.InflaT   = InflaVar(1);
    expParam.InflaV   = InflaVar(2);
    
    [mic_reSamp] = OfflineLoadBaselineVoice(dirs);
   
    DAQin = []; rawData = [];
    for ii = 1:expParam.numTrial
        expParam.curTrialNum  = ii;
        expParam.curTrial     = ['Trial' num2str(ii)];
        expParam.curSessTrial = [expParam.subject expParam.run expParam.curTrial];

        %Level of f0 change based on results from Laryngeal pert Exp
        audStimP = dfSetAudapFiles(dirs, expParam, ii);
        p.timeDomainPitchShiftSchedule = audStimP.pertSched;            

        %Cue to begin trial
        pause(expParam.cuePause)
        
        %Phonation Start
        fprintf('Trial %d\n', ii)
        AudapterIO('init', p);
        Audapter('reset');

        % Load the PreRecorded Baseline Mic signal
        % Split the signal into frames
        mic_frames = makecell(mic_reSamp, expParam.frameLen);

        for n = 1:length(mic_frames)
            Audapter('runFrame', mic_frames{n});
        end

        dataDAQ = [];

        %Save the data
        data    = AudapterIO('getData'); % This will need to become a try statement again
        DAQin   = cat(3, DAQin, dataDAQ);
        rawData = cat(1, rawData, data);
        
        dfSaveWavRec(data, expParam, dirs);
        
        %Play the sound, if you care to
        if expParam.bPlay; soundsc(data.signalOut, expParam.sRateAnal); end

        pause(10)
    end
    close all;
    elapsed_time = toc(ET)/60;    % Elapsed Time of the experiment
    fprintf('\nElapsed Time: %f (min)\n', elapsed_time)
    
    % Store all the variables and data from the session in a large structure
    expParam.elapsedTime = elapsed_time;

    OA.dirs        = dirs;
    OA.expParam    = expParam;
    OA.p           = p;
    OA.audStimP    = audStimP;
    OA.DAQin       = DAQin;
    OA.rawData     = rawData;

    dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.subject expParam.run '.mat']);
    save(dirs.RecFileDir, 'OA')
else
    dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.subject expParam.run '.mat']);
    load(dirs.RecFileDir)
end
close all

[~, auRes] = dfAnalysisAudapter(dirs, OA.expParam, OA.rawData);

drawAudRespMeanTrial(auRes, dirs.SavResultsDir)
pause(2)
drawAudRespIndivTrial(auRes, dirs.SavResultsDir)
pause(2)
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

function [mic_reSamp] = OfflineLoadBaselineVoice(dirs)
%Making an extra function because I am extra
trial = 1;

%Load previously recorded voice sample to perturb
fprintf('Loading Previously Recorded Data Set %s\n\n', dirs.BaseFile)
load(dirs.BaseFile, 'DRF')
baseData = DRF.rawData(trial);

fs       = DRF.expParam.sRateAnal;
mic      = baseData.signalIn;
micEnd   = mic(48000:end);
mic      = [mic; micEnd];
downFact = baseData.params.downFact;
sr       = baseData.params.sr;

%Resample at 48000Hz
mic_reSamp = resample(mic, sr*downFact, fs);
mic_reSamp = mic_reSamp - mean(mic_reSamp);
end

function micproc = OfflineLoadBaselineVoiceWav(dirs, expParam)

inputWav     = 'W:\Experiments\Dissociating-Role-of-Feedback-in-Voice-Motor-Control\Pilot0\Aud1\wavFiles\Pilot0Aud1Trial1_micIn.wav';
[mic, micfs] = audioread(inputWav);

fs = expParam.sRate;

mic_resamp = resample(mic, fs, micfs);
micproc    = mic_resamp - mean(mic_resamp);
end

function visPSSsigs(data)

frameDur = data.params.frameLen / data.params.sr;
tAxis = 0 : frameDur : frameDur * (size(data.fmts ,1) - 1);
% get the formant plot bounds

figure('Position', [200, 200, 800, 400]);
[s, f, t]    = spectrogram(data.signalIn, 64, 48, 1024, data.params.sr);
[s2, f2, t2] = spectrogram(data.signalOut, 64, 48, 1024, data.params.sr);

subplot(211);
imagesc(t, f, 10 * log10(abs(s))); 
axis xy;
colormap jet;
hold on;
plot(tAxis, data.fmts(:, 1 : 2), 'w', 'LineWidth', 2);
set(gca, 'YLim', [0, 4000]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Original');

subplot(212);
imagesc(t2, f2, 10 * log10(abs(s2)));
axis xy;
colormap jet;
hold on;
plot(tAxis, data.fmts(:, 1 : 2), 'w', 'LineWidth', 2);
plot(tAxis, data.sfmts(:, 1 : 2), 'g--', 'LineWidth', 2);
set(gca, 'YLim', [0, 4000]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Shifted');

figure;
hold on;
sigTAxis1 = 0 : 1 / data.params.sr : ...
    1 / data.params.sr * (length(data.signalIn) - 1);
plot(sigTAxis1, data.signalIn, 'b-');
ylabel('Input waveform');
xlabel('Time (s)');

sigTAxis2 = 0 : 1 / data.params.sr : ...
    1 / data.params.sr * (length(data.signalOut) - 1);
plot(sigTAxis2, data.signalOut, 'r-');
legend({'Input', 'Output'});
ylabel('Output waveform');
xlabel('Time (s)');

figure;
hold on;
plot(sigTAxis1(1 : end - 1), diff(data.signalIn), 'b-');
ylabel('d Input waveform');
xlabel('Time (s)');
plot(sigTAxis2(1 : end - 1), diff(data.signalOut), 'r-');
legend({'Input', 'Output'});
ylabel('d Output waveform');
xlabel('Time (s)');

% When time-domain pitch shifting is activated, data.pitchHz and
% data.shfitedPitchHz record the tracked and shifted pitch values in Hz,
% for supra-threshold frames.
figure;
plot(tAxis, data.pitchHz, 'b-');
hold on;
plot(tAxis, data.shiftedPitchHz, 'r-');
legend({'Input', 'Output'});
xlabel('Time (s)');
ylabel('Pitch (Hz)');

% axis([0 4 190 240])
end

function [f0b, targRMS, rmsB, gender] = loadBaselineVoice(dirs)

if exist(dirs.BaseFile, 'File')
    load(dirs.BaseFile, 'DRF')
    
    f0b     = DRF.qRes.meanf0;
    targRMS = DRF.qRes.meanRMS;
    rmsB    = DRF.expParam.rmsB;
    gender  = DRF.expParam.gender;
else
    fprintf('Could not find baseline voice file at %s\n', dirs.BaseFile)
    fprintf('Loading Default Values for f0b, meanRMS, and rmsB\n')
    f0b     = 100;
    targRMS = 70.00;
    rmsB    = 0.00002;
    gender  = 'female';
end
end