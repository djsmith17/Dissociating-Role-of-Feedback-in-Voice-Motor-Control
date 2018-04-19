function dfRunAFPerturbOffline()
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
debug = 0;

% Main Experimental prompt: Subject/Run Information
subject    = 'Pilot0';    % Subject#, Pilot#, null
run        = 'AF1';     % AF1, DS1, etc
blLoudness = 79.34;     % (dB SPL) Baseline loudness
gender     = 'male';  % "male" or "female"
InflaVarNm = 'IV1';
BaseRun    = 'BV1';
collectNewData         = 1; %Boolean

% Dialogue box asking for what type of Pitch-Shifted Feedback?
pertType = questdlg('What type of Perturbation?', 'Type of Perturbation?', 'Linear Standard', 'Sigmoid Matched', 'Sigmoid Matched');
switch pertType
    case 'Linear Standard'
        pertTypeSw = 0;
    case 'Sigmoid Matched'
        pertTypeSw = 1;
end

% Dialogue box asking if Practice set or Full set of trials
recType = questdlg('Practice or Full?','Length', 'Practice', 'Diagnostic', 'Full','Full');
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

%Experiment Configurations
expParam.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType      = 'Auditory Perturbation_Perceptual';
expParam.subject      = subject;
expParam.run          = [run 'Offline'];
expParam.curSess      = [expParam.subject expParam.run];
expParam.targRMS      = blLoudness;
expParam.gender       = gender;
expParam.balloon      = 'N/A';
expParam.tightness    = 'N/A';
expParam.InflaVarNm   = InflaVarNm;
expParam.niDev        = 'Dev2';          % NIDAQ Device Name. For more information, see dfInitNIDAQ
expParam.trialLen     = 4;               % Seconds
expParam.numTrial     = numTrials;
expParam.curTrial     = [];
expParam.perCatch     = perCatch;
expParam.AudFB        = 'Voice Shifted';
expParam.AudFBSw      = 1; %Voice Shifted
expParam.AudPert      = pertType;
expParam.AudPertSw    = pertTypeSw;
expParam.bVis         = 1;
expParam.bPlay        = 0;

expParam.baseRun      = BaseRun;
expParam.baseFile     = [expParam.subject expParam.baseRun 'DRF.mat'];

%Set our dirs based on the project
dirs = dfDirs(expParam.project);
% Folder paths to save data files
dirs.RecFileDir = fullfile(dirs.RecData, expParam.subject, expParam.run);
dirs.RecWaveDir = fullfile(dirs.RecFileDir, 'wavFiles');

if exist(dirs.RecFileDir, 'dir') == 0
    mkdir(dirs.RecFileDir)
end
if exist(dirs.RecWaveDir, 'dir') == 0
    mkdir(dirs.RecWaveDir)
end

% Look for the Inflation Response Files
expParam.InflaFile    = [expParam.subject expParam.InflaVarNm 'DRF.mat'];
dirs.InflaVarFile = fullfile(dirs.SavData, expParam.subject, expParam.InflaVarNm, expParam.InflaFile);
if ~exist(dirs.InflaVarFile, 'file')
    fprintf('ERROR: No Inflation Vars File at %s!\n', dirs.InflaVarFile)
    return
else
    fprintf('Inflation Variables found!!\n')
end

% Look for the Baseline Wav Files
dirs.SavBaseFile = fullfile(dirs.SavData, expParam.subject, expParam.baseRun, expParam.baseFile);
if ~exist(dirs.SavBaseFile, 'file')
    fprintf('ERROR: No voice file at %s!\n', dirs.SavBaseFile)
    return
end

% Look for a place to save the data
dirs.SavResultsDir = fullfile(dirs.Results, expParam.subject, expParam.run);
if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

if collectNewData == 1
    %Paradigm Configurations
    expParam.sRate              = 48000;  % Hardware sampling rate (before downsampling)
    expParam.frameLen           = 96;  % Before downsampling
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
    expParam.ostFN = fullfile(dirs.Prelim, 'AFPerturbOST.ost'); check_file(expParam.ostFN);
    expParam.pcfFN = fullfile(dirs.Prelim, 'AFPerturbPCF.pcf'); check_file(expParam.pcfFN);
    
    %Set up Auditory Feedback (Masking Noise, Pitch-Shift?)
    [expParam, p]      = dfSetAudFB(expParam, dirs, p);
    
    %Set up the order of trials (Order of perturbed, control, etc)
    expParam.trialType = dfSetTrialOrder(expParam.numTrial, expParam.perCatch); %numTrials, percentCatch

    %Select the trigger points for perturbation onset and offset and creating
    %the digital signal to be sent to the NIDAQ
    [expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType, 1);

    expParam.cuePause  = 1.0; % How long the cue period lasts
    expParam.buffPause = 0.2; % Give them a moment to start speaking
    expParam.resPause  = 2.0; % How long the rest/VisFB lasts
    expParam.boundsRMS = 3;   % +/- dB
    
    % Gives variable of InflaVar. Analyzed from previous recording
    load(dirs.InflaVarFile);
    expParam.InflaT   = InflaVar(1);
    expParam.InflaV   = InflaVar(2);
   
    % Load the PreRecorded Baseline Mic signal
    [mic_frames, f0b] = OfflineLoadBaselineVoice(dirs);
    expParam.f0b = f0b;

    DAQin = []; rawData = [];
    for ii = 1:expParam.numTrial
        expParam.curTrialNum  = ii;
        expParam.curTrial     = ['Trial' num2str(ii)];
        expParam.curSessTrial = [expParam.subject expParam.run expParam.curTrial];

        %Level of f0 change based on results from Laryngeal pert Exp
        audStimP = dfSetAudapFiles(expParam, dirs, ii, debug);

        %Set the OST and PCF functions
        Audapter('ost', expParam.ostFN, 0);
        Audapter('pcf', expParam.pcfFN, 0);
        
        %Setup which perturb file we want
        NIDAQsig = [expParam.sigs(:,ii) nVS];
        queueOutputData(s, NIDAQsig);
        
        %Cue to begin trial
        pause(expParam.cuePause)
        
        %Phonation Start

        fprintf('Trial %d\n', ii)
        AudapterIO('init', p);
        Audapter('reset');
        pause(expParam.buffPause)

        for n = 1:length(mic_frames)
            Audapter('runFrame', mic_frames{n});
        end

        [dataDAQ, ~] = s.startForeground;

        %Save the data
        data    = AudapterIO('getData'); % This will need to become a try statement again
        DAQin   = cat(3, DAQin, dataDAQ);
        rawData = cat(1, rawData, data);
        
        switch recType
            case 'Diagnostic'
                dfSaveWavRec(data, expParam, dirs);
            case 'Full'
                dfSaveWavRec(data, expParam, dirs);
        end
        
        %Play the sound, if you care to
        if expParam.bPlay; soundsc(data.signalIn, data.expParam.sRateAnal); end

        pause(expParam.resPause)
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

f0b = OA.expParam.f0b;
aFa = 1; iRf = 0;
niAn = struct;
niAn.sRate = 8000;
[auAn, auRes] = dfAnalysisAudapter(dirs, OA.expParam, OA.rawData, f0b, aFa, iRf, niAn);

drawAudRespMeanTrial(auRes, dirs.SavResultsDir)
drawAudRespIndivTrial(auRes, dirs.SavResultsDir)
end

function [mic_frames, f0b] = OfflineLoadBaselineVoice(dirs)
%Making an extra function because I am extra
trial = 1;

%Load previously recorded voice sample to perturb
fprintf('Loading Previously Recorded Data Set %s\n\n', dirs.SavBaseFile)
load(dirs.SavBaseFile);
baseData = DRF.rawData(trial);

fs       = DRF.expParam.sRateAnal;
mic      = [baseData.signalIn; zeros(fs*0.15,1)];
downFact = baseData.params.downFact;
sr       = baseData.params.sr;
frameLen = baseData.params.frameLen;

%Resample at 48000Hz
mic_reSamp = resample(mic, sr*downFact, fs);
%Split the signal into frames
mic_frames = makecell(mic_reSamp, frameLen*downFact);

f0b = DRF.qRes.meanf0;
end