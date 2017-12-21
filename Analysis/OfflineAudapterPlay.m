function OfflineAudapterPlay(varargin)
%This scripts loads a previously recorded audio signal and provides a
%Pitch-shift to it in a similiar fashion that happens during online
%testing.

%This script calls the following (8) functions:
%dfDirs.m
%initNIDAQ.m
%dfSetAudFB.m
%dfSetTrialOrder.m
%dfMakePertSignal.m
%dfSetAudapFiles.m

close all;
ET = tic;
rng('shuffle');

prompt = {'Subject ID:',...
          'Session ID:',...
          'Baseline Loudness (dB SPL):',...
          'Gender ("male" or "female"):',...
          'Inflation Vars:',...
          'Baseline Run:'};
name = 'Subject Information';
numlines = 1;
defaultanswer = {'null', 'AF1', '60', 'female', 'IV1', 'BV1'};
answer = inputdlg(prompt, name, numlines, defaultanswer);

if isempty(answer)
    return
end

pertType = questdlg('What type of Perturbation?', 'Type of Perturbation?', 'Linear Standard', 'Sinusoid Matched', 'Sinusoid Matched');
switch pertType
    case 'Linear Standard'
        pertTypeSw = 0;
    case 'Sinusoid Matched'
        pertTypeSw = 1;
end

collectNewData         = 0; %Boolean

%Experiment Configurations
expParam.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType      = 'Auditory Perturbation_Perceptual';
expParam.subject      = answer{1};
expParam.run          = [answer{2} 'Offline'];
expParam.targRMS      = str2double(answer{3});
expParam.gender       = answer{4};
expParam.curSess      = [expParam.subject expParam.run];
expParam.InflaVar     = answer{5};
expParam.baseRun      = answer{6};
expParam.numTrial     = 1;
expParam.curTrial     = [];
expParam.perCatch     = 1.00;
expParam.AudFB        = 'Voice Shifted';
expParam.AudFBSw      = 1; %Voice Shifted
expParam.trialLen     = 4; %Seconds
expParam.niDev        = 'Dev2';
expParam.bVis         = 1;
expParam.bPlay        = 0;
expParam.AudPert      = pertType;
expParam.AudPertSw    = pertTypeSw;

expParam.InflaFile    = [expParam.subject expParam.InflaVar 'DRF.mat'];
expParam.baseFile     = [expParam.subject expParam.baseRun 'DRF.mat'];

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

dirs.InflaVarFile = fullfile(dirs.SavData, expParam.subject, expParam.InflaVar, expParam.InflaFile);
if ~exist(dirs.InflaVarFile, 'file')
    fprintf('ERROR: No Inflation Vars File at %s!\n', dirs.InflaVarFile)
    return
end

dirs.SavBaseFile = fullfile(dirs.SavData, expParam.subject, expParam.baseRun, expParam.baseFile);
if ~exist(dirs.SavBaseFile, 'file')
    fprintf('ERROR: No voice file at %s!\n', dirs.SavBasFile)
    return
end

dirs.SavResultsDir = fullfile(dirs.Results, expParam.subject, expParam.run);
if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

if collectNewData == 1
    %Paradigm Configurations
    expParam.sRate              = 48000;  % Hardware sampling rate (before downsampling)
    expParam.downFact           = 3;
    expParam.sRateAnal          = expParam.sRate/expParam.downFact; %Everything get automatically downsampled! So annoying
    expParam.frameLen           = 96;  % Before downsampling
    expParam.audioInterfaceName = 'MOTU MicroBook'; %'ASIO4ALL' 'Komplete'

    %Set up Audapter
    Audapter('deviceName', expParam.audioInterfaceName);
    Audapter('setParam', 'downFact', expParam.downFact, 0);
    Audapter('setParam', 'sRate', expParam.sRateAnal, 0);
    Audapter('setParam', 'frameLen', expParam.frameLen / expParam.downFact, 0);
    p = getAudapterDefaultParams(expParam.gender);

    %Set up Parameters to control NIDAQ and Perturbatron
    [s, niCh, nVS]  = initNIDAQ(expParam.niDev, expParam.trialLen);
    expParam.sRateQ = s.Rate; % NIDAQ sampling rate
    expParam.niCh   = niCh;   % Structure of Channel Names

    %Set up OST and PCF Files
    expParam.ostFN = fullfile(dirs.Prelim, 'AFPerturbOST.ost'); check_file(expParam.ostFN);
    expParam.pcfFN = fullfile(dirs.Prelim, 'AFPerturbPCF.pcf'); check_file(expParam.pcfFN);
    
    [expParam, p]      = dfSetAudFB(expParam, dirs, p); %Trials with masking or no...
    
    expParam.trialType = dfSetTrialOrder(expParam.numTrial, expParam.perCatch); %numTrials, percentCatch

    [expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType, 1);

    expParam.cuePause  = 1.0;
    expParam.resPause  = 2.0;
    expParam.boundsRMS = 3;  %+/- dB
    
    % Gives variable of InflaVar. Analyzed from previous recording
    load(dirs.InflaVarFile);
    expParam.InflaT   = InflaVar(1);
    expParam.InflaV   = InflaVar(2);
   
    % Load the PreRecorded Baseline Mic signal
    [mic_frames] = OfflineLoadBaselineVoice(dirs);

    DAQin = []; rawData = [];
    for ii = 1:expParam.numTrial
        expParam.curTrial     = ['Trial' num2str(ii)];
        expParam.curSessTrial = [expParam.subject expParam.run expParam.curTrial];

        %Level of f0 change based on results from Laryngeal pert Exp
        audStimP = dfSetAudapFiles(expParam, dirs, ii, 0);

        %Set the OST and PCF functions
        Audapter('ost', expParam.ostFN, 0);
        Audapter('pcf', expParam.pcfFN, 0);
        
        %Setup which perturb file we want
        NIDAQsig = [expParam.sigs(:,ii) nVS];
        queueOutputData(s, NIDAQsig);

        fprintf('Trial %d\n', ii)
        AudapterIO('init', p);
        Audapter('reset');

        for n = 1:length(mic_frames)
            Audapter('runFrame', mic_frames{n});
        end

        [dataDAQ, time] = s.startForeground;

        %Save the data
        data    = dfSaveRawData(expParam, dirs);
        DAQin   = cat(3, DAQin, dataDAQ);
        rawData = cat(1, rawData, data);
        
        %Play the sound, if you care to
        if expParam.bPlay; soundsc(data.signalIn, data.expParam.sRateAnal); end

        pause(expParam.resPause)
    end
    
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

GTFile = [expParam.subject 'GT1' 'DRF.mat'];
dirs.GTBaseFile = fullfile(dirs.SavData, expParam.subject, 'GT1', GTFile);
if ~exist(dirs.SavBaseFile, 'file')
    fprintf('ERROR: No GT file at %s!\n', dirs.GTBasFile)
    return
end
load(dirs.GTBaseFile)

bTf0b = GT.subjf0;
[auAn, auRes] = dfAnalysisAudapter(dirs, OA.expParam, OA.rawData, bTf0b, 1);

drawAudRespMeanTrial(auRes, dirs.SavResultsDir)
end

function [mic_frames] = OfflineLoadBaselineVoice(dirs)
%Making an extra function because I am extra
trial = 2;

%Load previously recorded voice sample to perturb
fprintf('Loading Previously Recorded Data Set %s\n\n', dirs.SavBaseFile)
load(dirs.SavBaseFile);
baseData = DRF.rawData(trial);

mic      = baseData.signalIn;
fs       = DRF.expParam.sRateAnal;
downFact = baseData.params.downFact;
sr       = baseData.params.sr;
frameLen = baseData.params.frameLen;

%Resample at 48000Hz
mic_reSamp = resample(mic, sr*downFact, fs);
%Split the signal into frames
mic_frames = makecell(mic_reSamp, frameLen*downFact);
end