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
          'Gender ("male" or "female"):'};
name = 'Subject Information';
numlines = 1;
defaultanswer = {'null', 'AF1', '60', 'female'};
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
expParam.run          = answer{2};
expParam.targRMS      = str2double(answer{3});
expParam.gender       = answer{4};
expParam.curSess      = [expParam.subject expParam.run];
expParam.numTrial     = 10;
expParam.curTrial     = [];
expParam.perCatch     = 1.00;
expParam.AudFB        = 'Voice Shifted';
expParam.AudFBSw      = 1; %Voice Shifted
expParam.trialLen     = 4; %Seconds
expParam.niDev        = 'Dev2';
expParam.bVis         = 0;
expParam.AudPert      = pertType;
expParam.AudPertSw    = pertTypeSw;

expParam.baseRec      = 'BV1';

expParam.bf0Vis        = 0;
expParam.bPlay         = 0;
expParam.offLineTrial  = 37;

dirs = dfDirs(expParam.project);
% Folder paths to save data files
dirs.RecFileDir  = fullfile(dirs.RecData, expParam.subject, [expParam.run 'Offline']);
dirs.RecWaveDir  = fullfile(dirs.RecFileDir, 'wavFiles');

dirs.SavFileDir    = fullfile(dirs.SavData, expParam.subject, expParam.baseRec);
dirs.SavResultsDir = fullfile(dirs.Results, expParam.subject, 'offline');
dirs.saveFileSuffix = '_offlinePSR';

if exist(dirs.RecFileDir, 'dir') == 0
    mkdir(dirs.RecFileDir)
end
if exist(dirs.RecWaveDir, 'dir') == 0
    mkdir(dirs.RecWaveDir)
end
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
    s = initNIDAQ(expParam.trialLen, 'Dev2');
    expParam.sRateQ = s.Rate; %save the sampling rate of the NIDAQ

    %Set up OST and PCF Files
    expParam.ostFN = fullfile(dirs.Prelim, 'AFPerturbOST.ost'); check_file(expParam.ostFN);
    expParam.pcfFN = fullfile(dirs.Prelim, 'AFPerturbPCF.pcf'); check_file(expParam.pcfFN);

    %Should return variables of InflaRespRoute and tStep. 
    %Recorded from previous experiments
    dirs.InflaRespFile = fullfile(dirs.SavData, expParam.subject, [expParam.subject '_AveInflaResp.mat']);
    try
        load(dirs.InflaRespFile);
        expParam.InflaRespRoute = InflaRespRoute;
        expParam.tStep          = tStep;
    catch me
        fprintf('\nSubject Data does not exist at %s \n', dirs.InflaRespFile)
    end

    [expParam, p]      = dfSetAudFB(expParam, dirs, p); %Trials with masking or no... 

    expParam.trialType = dfSetTrialOrder(expParam.numTrial, expParam.perCatch); %numTrials, percentCatch

    [expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType);

    %Create a negative voltage signal for the force sensors
    negVolSrc = zeros(expParam.sRateQ*expParam.trialLen, 1) - 1;
    negVolSrc(1) = 0; negVolSrc(end) = 0;

    expParam.resPause = 2.0;

    %Taking the first trial for ease. File out will be 'data'
    fprintf('Loading Previously Recorded Data Set...\n\n')
    d = dir([dirs.SavFileDir, '\*.mat']);
    fnames = sort_nat({d.name}); 
    load(fullfile(dirs.SavFileDir, fnames{1})); 
    data = DRF.rawData(expParam.offLineTrial);
    
    Mraw  = data.signalIn; 
    fs    = data.params.sRate;

    DAQin   = [];
    rawData = [];
    for ii = 1:expParam.numTrial
        expParam.curTrial    = ['Trial' num2str(ii)];
        expParam.curExpTrial = [expParam.subject expParam.run expParam.curTrial];

        audStimP = dfSetAudapFiles(expParam, dirs, ii, 1);

        %Set the OST and PCF functions
        Audapter('ost', expParam.ostFN, 0);
        Audapter('pcf', expParam.pcfFN, 0);

        %Resample at 48000Hz
        Mraw_reSamp = resample(Mraw, data.params.sr * data.params.downFact, fs);
        %Split the signal into frames
        Mraw_frames = makecell(Mraw_reSamp, data.params.frameLen * data.params.downFact);

        fprintf('Running Trial %d\n', ii)
        AudapterIO('init', p);
        Audapter('reset');

        for n = 1:length(Mraw_frames)
            Audapter('runFrame', Mraw_frames{n});
        end
        
        NIDAQsig = [expParam.sigs(:,ii) negVolSrc];
        queueOutputData(s, NIDAQsig);

        [dataDAQ, time] = s.startForeground;

        data = dfSaveRawData(expParam, dirs);
        DAQin = cat(3, DAQin, dataDAQ);
        rawData = cat(1, rawData, data);
        
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

[auAn, res] = dfAnalysisAudapter(OA.expParam, OA.rawData, OA.DAQin);

drawAudResp_AllTrial(res, auAn.curSess, OA.expParam.curRec, dirs.SavResultsDir)

drawAudResp_InterTrial(res.timeSec, res.meanTrialf0_St, res.meanTrialf0_Sp, res.f0LimitsSec, res.trialCount, res.meanTrialf0b, auAn.curSess, OA.expParam.curRec, dirs.SavResultsDir)
end