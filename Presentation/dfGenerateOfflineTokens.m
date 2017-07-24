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

expParam.trialType = dfSetTrialOrder(expParam.numTrial, expParam.per); %numTrials, percentCatch
expParam.trigs = [0 expParam.trialLen];

% [expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType);

audStimP = dfSetAudapFiles(expParam.ostFN, expParam.pcfFN, expParam.trialType, expParam.trigs, expParam.stimType, InflaRespRoute, tStep);

%Set the OST and PCF functions
Audapter('ost', expParam.ostFN, 0);
Audapter('pcf', expParam.pcfFN, 0);

fprintf('Running Trial %d\n', ii)
AudapterIO('init', p);
Audapter('reset');

for n = 1:length(expParam.baseToken_frames)
    Audapter('runFrame', expParam.baseToken_frames{n});
end

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