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

pertType = questdlg('What type of Perturbation?', 'Type of Perturbation?', '-100 cents ramped', 'Laryngeal Pert Matched', 'Laryngeal Pert Matched');
switch pertType
    case '-100 cents ramped'
        pertTypeSw = 0;
    case 'Laryngeal Pert Matched'
        pertTypeSw = 1;
end

collectNewData         = 1; %Boolean

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
expParam.numTrial     = 10;
expParam.curTrial     = [];
expParam.perCatch     = 1.00;
expParam.AudFB        = 'Voice Shifted';
expParam.AudFBSw      = 1; %Voice Shifted
expParam.trialLen     = 4; %Seconds
expParam.niDev        = 'Dev2';
expParam.bVis         = 0;
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

dirs.SavResultsDir = fullfile(dirs.Results, expParam.subject, 'offline');
dirs.saveFileSuffix = '_offlinePSR';

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

    [expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType);

    expParam.cuePause  = 1.0;
    expParam.resPause  = 2.0;
    expParam.boundsRMS = 3;  %+/- dB
    
    % Gives variable of InflaVar. Analyzed from previous recording
    load(dirs.InflaVarFile);
    expParam.InflaT   = InflaVar(1);
    expParam.InflaV   = InflaVar(2);
   
    % Load the PreRecorded Baseline Mic signal
    [PreRMic, PreRfs] = OfflineLoadBaselineVoice(dirs);

    DAQin = []; rawData = [];
    for ii = 1:expParam.numTrial
        expParam.curTrial    = ['Trial' num2str(ii)];
        expParam.curExpTrial = [expParam.subject expParam.run expParam.curTrial];

        %Level of f0 change based on results from Laryngeal pert Exp
        audStimP = dfSetAudapFiles(expParam, dirs, ii, 1);

        %Set the OST and PCF functions
        Audapter('ost', expParam.ostFN, 0);
        Audapter('pcf', expParam.pcfFN, 0);
        
        %Setup which perturb file we want
        NIDAQsig = [expParam.sigs(:,ii) nVS];
        queueOutputData(s, NIDAQsig);

        %Resample at 48000Hz
        mic_reSamp = resample(PreRMic, data.params.sr * data.params.downFact, PreRfs);
        %Split the signal into frames
        mic_frames = makecell(mic_reSamp, data.params.frameLen * data.params.downFact);

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
        if expParam.bPlay; soundsc(data_off.signalIn, data_off.expParam.sRateAnal); end

        pause(expParam.resPause)
    end
    
    OA.dirs        = dirs;
    OA.expParam    = expParam;
    OA.p           = p;
    OA.audStimP    = audStimP;
    OA.DAQin       = DAQin;
    OA.rawData     = rawData;      

    dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.subject dirs.saveFileSuffix '.mat']);
    save(dirs.RecFileDir, 'OA')
else
    dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.subject dirs.saveFileSuffix '.mat']);
    load(dirs.RecFileDir)
end
close all

% [auAn, res] = dfAnalysisAudapter(OA.expParam, OA.rawData, OA.DAQin);
% 
% drawAudResp_AllTrial(res, auAn.curSess, OA.expParam.curRec, dirs.SavResultsDir)
% 
% drawAudResp_InterTrial(res.timeSec, res.meanTrialf0_St, res.meanTrialf0_Sp, res.f0LimitsSec, res.trialCount, res.meanTrialf0b, auAn.curSess, OA.expParam.curRec, dirs.SavResultsDir)
end

function [mic, fs] = OfflineLoadBaselineVoice(dirs)
%Making an extra function because I am extra
trial = 2;

%Load previously recorded voice sample to perturb
fprintf('Loading Previously Recorded Data Set %s\n\n', dirs.SavBaseFile)
load(dirs.SavBaseFile);
baseData = DRF.rawData(trial);

mic = baseData.signalIn;
fs  = DRF.expParam.sRateAnal;
end