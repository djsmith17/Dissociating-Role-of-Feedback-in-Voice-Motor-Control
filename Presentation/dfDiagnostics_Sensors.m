function dfDiagnostics_Sensors(varargin)
%A quick test of the force sensors before running the actual experiment.
%This makes sure that the sensors are working they should be and we can
%continue with the experiment. Eventually this will also include the
%pressure sensor. 

%This script calls the following (4) functions:
%dfDirs.m
%initNIDAQ.m
%dfMakePertSignal.m
%drawDAQsignal.m
close all;

prompt = {'Subject ID:',...
          'Session ID:',...
          'Number of Trials:',...
          'Percent Perturbed (Dec)',...
          'Collect New Data?:'};
name = 'Subject Information';
numlines = 1;
defaultanswer = {'null', 'DS1', '5', '1', 'yes'};
answer = inputdlg(prompt, name, numlines, defaultanswer);

if isempty(answer)
    return
end

expParam.project       = 'NIDAQSensorDiagnostics';
expParam.expType       = 'Somatosensory Perturbation_Perceptual';
expParam.subject       = answer{1}; %Subject#, Pilot#, null
expParam.run           = answer{2};
expParam.curSess       = [expParam.subject ' ' expParam.run];
expParam.numTrial      = str2double(answer{3});
expParam.perCatch      = str2double(answer{4});
expParam.trialLen      = 4; %Seconds
expParam.niDev         = 'Dev2';
expParam.AudFB         = 'Masking Noise';
expParam.gender        = 'N/A';
sv2F                   = 1; %Boolean
collectNewData         = answer{5};

dirs = dfDirs(expParam.project);

dirs.RecFileDir    = fullfile(dirs.RecData, expParam.subject, expParam.run);
dirs.SavResultsDir = fullfile(dirs.RecData, expParam.subject, expParam.run); %Where to save results 

if exist(dirs.RecFileDir, 'dir') == 0
    mkdir(dirs.RecFileDir)
end
if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end
dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.subject expParam.run 'NSD.mat']);

if strcmp(collectNewData, 'yes')
    expParam.sRate       = 48000;
    expParam.downFact    = 3;
    expParam.sRateAnal   = expParam.sRate/expParam.downFact; %Everything get automatically downsampled! So annoying
    expParam.resPause    = 3;
    
    [s, niCh, nVS]  = initNIDAQ(expParam.niDev, expParam.trialLen);
    expParam.sRateQ = s.Rate;
    expParam.niCh   = niCh;

    expParam.trialType              = dfSetTrialOrder(expParam.numTrial, expParam.perCatch);
    [expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType, 1);  

    DAQin = []; DAQtime = [];
    for ii = 1:expParam.numTrial
        NIDAQsig = [expParam.sigs(:,ii) nVS];
        queueOutputData(s, NIDAQsig);
        fprintf('Running Trial %d\n', ii)
        
        tic
        [data_DAQ, time] = s.startForeground;
        toc

        DAQin   = cat(3, DAQin, data_DAQ);
        DAQtime = cat(3, DAQtime, time);
        
        pause(expParam.resPause)      
    end
    
    NSD.expParam    = expParam;
    NSD.dirs        = dirs;
    NSD.DAQin       = DAQin;

    save(dirs.RecFileDir, 'NSD')
else
    load(dirs.RecFileDir)
end

f0b = 100;
pF  = 1;
iRF = 0;
[niAn, niRes] = dfAnalysisNIDAQ(dirs, NSD.expParam, NSD.DAQin, f0b, 0, iRF, pF);

niRes.numPertTrialsNi = niRes.numPertTrials;

% drawDAQsignal(niAn, 1, dirs.SavResultsDir, sv2F)
drawDAQAlignedPressure(niRes, dirs.SavResultsDir, sv2F)
% drawDAQAll(niAn, dirs.SavResultsDir, sv2F)
% drawDAQPresMic(niAn, dirs.SavResultsDir, sv2F)
end

function res = combineRes(niRes, auRes)

res.expType = auRes.expType;
res.subject = auRes.subject;
res.run     = auRes.run;
res.curSess = auRes.curSess;
res.AudFB   = auRes.AudFB;

res.f0Type  = auRes.f0Type;
res.etMH    = auRes.etMH;

res.numTrials     = auRes.numTrials;
res.audioM        = auRes.audioM;
res.audioH        = auRes.audioH;
res.svIdx         = auRes.svIdx;
res.expTrigsSv    = auRes.expTrigsSv;
res.pertIdx       = auRes.pertIdx;     % The indices of the svIdx;
res.pertTrig      = auRes.pertTrig;
res.contIdx       = auRes.contIdx;     % The indices of the svIdx;
res.contTrig      = auRes.contTrig;
res.numSaveTrials = auRes.numSaveTrials;
res.numContTrials = auRes.numContTrials;
res.numPertTrials = auRes.numPertTrials;

res.numPertTrialsNi = niRes.numPertTrials;
res.numContTrialsNi = niRes.numContTrials;
res.timeS           = niRes.timeS;
res.sensorP         = niRes.sensorP;        % Individual Processed perturbed trials. 
res.lagTimeP        = niRes.lagTimeP;
res.lagTimePm       = niRes.lagTimePm;
res.riseTimeP       = niRes.riseTimeP;
res.riseTimePm      = niRes.riseTimePm;
res.OnOfValP        = niRes.OnOfValP;
res.OnOfValPm       = niRes.OnOfValPm;
res.limitsP         = niRes.limitsP;

res.timeSAl       = niRes.timeSAl;
res.sensorPAl     = niRes.sensorPAl;
res.limitsPAl     = niRes.limitsPAl;

res.timef0        = auRes.timef0;
res.f0b           = auRes.f0b;

res.numContTrialsPP = auRes.numContTrialsPP;
res.numPertTrialsPP = auRes.numPertTrialsPP;
res.pertTrigPP      = auRes.pertTrigPP;

%Full Individual Trials: Mic/Head f0 Trace 
res.audioMf0TrialPert = auRes.audioMf0TrialPert;
res.audioMf0TrialCont = auRes.audioMf0TrialCont;
res.audioHf0TrialPert = auRes.audioHf0TrialPert;
res.audioHf0TrialCont = auRes.audioHf0TrialCont;
res.limitsA           = auRes.limitsA;
res.limitsAudRes      = auRes.limitsAudRes;

%Sections Trials: Mic/Head f0
res.secTime          = auRes.secTime;
res.audioMf0SecPert  = auRes.audioMf0SecPert;
res.audioMf0SecCont  = auRes.audioMf0SecCont;
res.audioHf0SecPert  = auRes.audioHf0SecPert;
res.audioHf0SecCont  = auRes.audioHf0SecCont;

%Mean Sectioned Trials: Mic/Head f0 Trace 
res.audioMf0MeanPert = auRes.audioMf0MeanPert; % [MeanSigOn 90%CI MeanSigOff 90%CI]
res.audioMf0MeanCont = auRes.audioMf0MeanCont;
res.audioHf0MeanPert = auRes.audioHf0MeanPert;
res.audioHf0MeanCont = auRes.audioHf0MeanCont;
res.limitsAmean      = auRes.limitsAmean;
res.limitsAMH        = auRes.limitsAMH;      % Limits Audio Corrected for MicXHead

%Inflation Response
res.respVar      = auRes.respVar;
res.respVarM     = auRes.respVarM;
res.respVarSD    = auRes.respVarSD;
res.InflaStimVar = auRes.InflaStimVar;

%NIAu Delay
res.allAuNiDelays = auRes.allAuNiDelays;

%Check the Ni trials against the Au trials
presInd = [];
if res.numPertTrialsPP < res.numPertTrialsNi
    setPertTrials = res.svIdx(res.pertIdx);
    for ii = 1:length(niRes.pertIdx)
        ind = [];
        ind = find(setPertTrials == niRes.pertIdx(ii));
        if ~isempty(ind)
            presInd = cat(1, presInd, ii);
        end
    end
    
    res.sensorPsv  = res.sensorP(:, presInd);
    res.lagTimePsv = res.lagTimeP(presInd, :);
    res.riseTimePsv= res.riseTimeP(presInd);
    res.OnOfValPsv = res.OnOfValP(presInd, :);
else
    res.sensorPsv    = res.sensorP;
    res.lagTimePsv   = res.lagTimeP;
    res.riseTimePsv  = res.riseTimeP;
    res.OnOfValPsv   = res.OnOfValP;   
end
end