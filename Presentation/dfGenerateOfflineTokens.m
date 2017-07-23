function modifiedToken = dfGenerateOfflineTokens(baseToken, dirs, gender, level)
%This scripts loads a previously recorded audio signal and provides a
%Pitch-shift to it in a similiar fashion that happens during online
%testing.

%Paradigm Configurations
expParam.baseToken          = baseToken;
expParam.gender             = gender;
expParam.trialLen           = 2; %seconds
expParam.numTrial           = 1;
expParam.per                = 1;

expParam.sRate              = 48000;  % Hardware sampling rate (before downsampling)
expParam.downFact           = 3;
expParam.sRateAnal          = expParam.sRate/expParam.downFact; %Everything get automatically downsampled! So annoying
expParam.frameLen           = 96;  % Before downsampling
expParam.audioInterfaceName = 'MOTU MicroBook'; %'ASIO4ALL' 'Komplete'

%Resample at 48000Hz
expParam.baseToken_reSamp = resample(expParam.baseToken, 48000, 16000);
%Split the signal into frames
expParam.baseToken_frames = makecell(expParam.baseToken_reSamp, expParam.frameLen);

%Set up Audapter
Audapter('deviceName', expParam.audioInterfaceName);
Audapter('setParam', 'downFact', expParam.downFact, 0);
Audapter('setParam', 'sRate', expParam.sRateAnal, 0);
Audapter('setParam', 'frameLen', expParam.frameLen / expParam.downFact, 0);
p = getAudapterDefaultParams(expParam.gender);

% %Set up Parameters to control NIDAQ and Perturbatron
% s = initNIDAQ(expParam.trialLen, 'Dev2');
% expParam.sRateQ = s.Rate; %save the sampling rate of the NIDAQ

%Set up OST and PCF Files
expParam.ostFN = fullfile(dirs.Prelim, 'AFPerturbOST.ost'); check_file(expParam.ostFN);
expParam.pcfFN = fullfile(dirs.Prelim, 'AFPerturbPCF.pcf'); check_file(expParam.pcfFN);

% %Should return variables of InflaRespRoute and tStep. 
% %Recorded from previous experiments
% dirs.InflaRespFile = fullfile(dirs.SavData, expParam.subject, [expParam.subject '_AveInflaResp.mat']);
% try
%     load(dirs.InflaRespFile);
% catch me
%     fprintf('\nSubject Data does not exist at %s \n', dirs.InflaRespFile)
% end

% [expParam, p]      = dfSetAudFB(expParam, dirs, p); %Trials with masking or no... 

expParam.trialType = dfSetTrialOrder(expParam.numTrial, expParam.per); %numTrials, percentCatch
expParam.trigs = [0 expParam.trialLen];

% [expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType);

%Create a negative voltage signal for the force sensors
% negVolSrc = zeros(expParam.sRateQ*expParam.trialLen, 1) - 1;
% negVolSrc(1) = 0; negVolSrc(end) = 0;

% expParam.resPause = 2.0;

%Taking the first trial for ease. File out will be 'data'
% fprintf('Loading Previously Recorded Data Set...\n\n')
% d = dir([dirs.SavFileDir, '\*.mat']);
% fnames = sort_nat({d.name}); 
% load(fullfile(dirs.SavFileDir, fnames{1})); 
% data = DRF.rawData(expParam.offLineTrial);

% Mraw  = data.signalIn; 
% fs    = data.params.sRate;

% DAQin   = [];
% rawData = [];
% for ii = 1:expParam.numTrial

%     expParam.curTrial   = ['Trial' num2str(ii)];
%     expParam.curExpTrial = [expParam.subject expParam.run expParam.curTrial];

audStimP = dfSetAudapFiles(InflaRespRoute, tStep, expParam.ostFN, expParam.pcfFN, expParam.trialType, expParam.trigs, expParam.stimType);

%Set the OST and PCF functions
Audapter('ost', expParam.ostFN, 0);
Audapter('pcf', expParam.pcfFN, 0);

fprintf('Running Trial %d\n', ii)
AudapterIO('init', p);
Audapter('reset');

for n = 1:length(expParam.baseToken_frames)
    Audapter('runFrame', expParam.baseToken_frames{n});
end

% NIDAQsig = [expParam.sigs(:,ii) negVolSrc];
% queueOutputData(s, NIDAQsig);

% [dataDAQ, time] = s.startForeground;

try
    data = AudapterIO('getData');
    modifiedToken = data.signalOut;
catch
    sprintf('\nAudapter decided not to show up today')
    data = [];
    modifiedToken = [];
    return
end

% data = dfSaveRawData(expParam, dirs);
% DAQin = cat(3, DAQin, dataDAQ);
% rawData = cat(1, rawData, data);

% if expParam.bPlay; soundsc(data_off.signalIn, data_off.expParam.sRateAnal); end

% pause(expParam.resPause)
% end

% OA.dirs        = dirs;
% OA.expParam    = expParam;
% OA.p           = p;
% OA.audStimP    = audStimP;
% OA.DAQin       = DAQin;
% OA.rawData     = rawData;      
% 
% dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.subject dirs.saveFileSuffix '.mat']);
% save(dirs.RecFileDir, 'OA')
end