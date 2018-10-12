function dfRunSubjAnalysisEndo()
% dfRunSubjAnalysisEndo() is my main script for analyzing recorded audio files 
% and sensor information from experiments studying Voice Motor Control. 
% This function is set up to analyze multiple subject and runs in an 
% identical fashion, so that once you have a new data set, you can run it 
% all very quickly. 
% Most importantly, these analyses calculate change in f0 of a speaker's 
% voice as they complete somatosensory and auditory feedback perturbation 
% tasks. 
%
% This makes use of the following functions:
% dfAnalysisNIDAQ.m
% dfAnalysisAudapter.m
%
% Requires the Signal Processing Toolbox

close all
AVar.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
AVar.participants  = {'DRF_ENP3'};    % List of multiple participants.
AVar.numPart       = length(AVar.participants);
AVar.run           = 'SFL1';
AVar.trials        = 10;
AVar.baselineFile  = 'BV1';            % Baseline Voice information
AVar.debug         = 0;
AVar.cond          = {'All'};
AVar.numCond       = length(AVar.cond);

dirs               = dfDirs(AVar.project);
dirs.LoadData      = dirs.SavData;

for i = 1:AVar.numPart
    participant = AVar.participants{i};
    run         = AVar.run;
    dirs.baselineData  = fullfile(dirs.LoadData, participant, AVar.baselineFile, [participant AVar.baselineFile 'DRF.mat']); % Where to find data
      
    % Look for the baseline voice info for this participant, then load it
    if exist(dirs.baselineData, 'file') == 0
        fprintf('ERROR: Could not find baseline data set at %s\n', dirs.baselineData)
        return
    else
        fprintf('Loading baseline data set for %s %s\n', participant, AVar.baselineFile)
        load(dirs.baselineData) % Returns DRF
        bV = DRF;
    end
    
    % Baseline f0
    f0b = bV.qRes.meanf0;
    
    % Define where to load raw data and save analyzed results
    dirs.SavResultsDir = fullfile(dirs.Results, participant, run);                                % Where to save results

    % Make sure there is a place to save results
    if exist(dirs.SavResultsDir, 'dir') == 0
        mkdir(dirs.SavResultsDir)
    end
     
    ssVF  = createSortingStruct();
    ssMN  = createSortingStruct();
    ssAll = createSortingStruct();
    ssAll.SeqAudFB = {}; ssAll.SeqAudFBSw = [];
    for j = 1:AVar.trials
        dirs.SavFileDir    = fullfile(dirs.LoadData, participant, run, [participant run 'Trial' num2str(j) 'DRF.mat']);  % Where to find data        
        
        % Look for the recording sessoin raw data for this participant, then load it
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%\n')
        if exist(dirs.SavFileDir, 'file') == 0
            fprintf('ERROR: Could not find saved data set at %s\n', dirs.SavFileDir)
            return
        else
            fprintf('Loading saved data set for %s %s\n', participant, run)
            load(dirs.SavFileDir) % Returns DRF
        end
        
        if strcmp(DRF.expParam.AudFB, 'Voice Feedback')
            [DRFVF, ssVF] = ammendSortingStruc(DRF, ssVF, j);
        else        
            [DRFMN, ssMN] = ammendSortingStruc(DRF, ssMN, j);
        end
        [DRFAll, ssAll] = ammendSortingStruc(DRF, ssAll, j);
        ssAll.SeqAudFB   = cat(1, ssAll.SeqAudFB, DRF.expParam.AudFB);
        ssAll.SeqAudFBSw = cat(1, ssAll.SeqAudFBSw, DRF.expParam.AudFBSw);
    end
    
    for ii = 1:AVar.numCond
        ext = AVar.cond{ii};
        ss  = eval(['ss' ext]);
        DRF = eval(['DRF' ext]);
        
        DRF.expParam.curSess   = [DRF.expParam.curSess ext];
        DRF.expParam.numTrial  = length(ss.trialType);
        DRF.expParam.trialType = ss.trialType;
        DRF.expParam.sigs      = ss.sigs;
        DRF.expParam.trigs     = ss.trigs;
        
        if strcmp(ext, 'All')
            DRF.expParam.SeqAudFB   = ss.SeqAudFB;
            DRF.expParam.SeqAudFBSw = ss.SeqAudFBSw;
        end            
        
        DRF.rawData = ss.rawData;
        DRF.DAQin   = ss.DAQin;
        
        % Identify the type of experiment and decide what types of analyzes
        % we need. pF: Pressure Flag; iRF: Inflation Response Flag
        [DRF, pF, iRF] = preAnalysisCheck(DRF, f0b);
        aFn = 0; aFa = 1; %Audio Analysis Flag

        % Analysis on the NIDAQ raw data
        [niAn, niRes] = dfAnalysisNIDAQ(dirs, DRF.expParam, DRF.DAQin, f0b, aFn, iRF, pF);
        % Analysis on the Audapter raw data
        [auAn, auRes] = dfAnalysisAudapter(dirs, DRF.expParam, DRF.rawData, f0b, aFa, iRF, niAn);

        % Combine Audapter and NIDAQ results into one neat MATLAB structure
        res = combineRes(niRes, auRes);

        dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [participant run ext 'ResultsDRF.mat']);
        if AVar.debug == 0
            % Save the results of this recording session
            fprintf('\nSaving Results for %s %s\n', participant, run)
            save(dirs.SavResultsFile, 'res')
        end
    end
        
end
end

function res = combineRes(niRes, auRes)

res.subject = auRes.subject;
res.run     = auRes.run;
res.curSess = auRes.curSess;
res.gender  = auRes.gender;
res.age     = auRes.age;

res.expType = auRes.expType;
res.AudFB   = auRes.AudFB;
res.SeqAudFB = auRes.SeqAudFB;

res.f0Type  = auRes.f0Type;
res.etMH    = auRes.etMH;

% Raw Recorded data
res.sRate               = auRes.sRate;
res.numTrial            = auRes.numTrial;   % Total trials recorded
res.removedTrialTracker = auRes.removedTrialTracker;
res.incTrialInfo        = auRes.incTrialInfo;
res.allAuMHDelays       = auRes.allAuMHDelays; % Vector of the delays between the NIDAQ and Audapter microphone recordings
res.allAuNiDelays       = auRes.allAuNiDelays; % Vector of the delays between the NIDAQ and Audapter microphone recordings

res.balloon         = niRes.balloon;
res.numPertTrialsNi = niRes.numPertTrials;
res.numContTrialsNi = niRes.numContTrials;
res.pertIdxNi       = niRes.pertIdx;

res.timeS           = niRes.timeS;
res.sensorP         = niRes.sensorP;        % Individual Processed perturbed trials. 
res.lagTimeP        = niRes.lagTimeP;
res.lagTimePm       = niRes.lagTimePm;
res.riseTimeP       = niRes.riseTimeP;
res.riseTimePm      = niRes.riseTimePm;
res.riseTimePSE     = niRes.riseTimePSE;
res.OnOfValP        = niRes.OnOfValP;
res.OnOfValPm       = niRes.OnOfValPm;
res.OnOfValPSE      = niRes.OnOfValPSE;
res.pTrialLossP     = niRes.pTrialLossP;
res.pTrialLossPm    = niRes.pTrialLossPm; 
res.pTrialLossPSE   = niRes.pTrialLossPSE;
res.limitsP         = niRes.limitsP;

% Sectioned and Aligned Pressure recordings 
res.timeSAl       = niRes.timeSAl;
res.sensorPAl     = niRes.sensorPAl;  
res.limitsPAl     = niRes.limitsPAl;

res.secTimeP    = niRes.secTimeP;
res.sensorPSec  = niRes.sensorPSec;
res.sensorPMean = niRes.sensorPMean;
res.limitsPMean = niRes.limitsPMean;

% Audio f0 analysis
res.timef0        = auRes.timef0;
res.f0b           = auRes.f0b;

% Which trials did I end up saving and plotting at the end of the day?
res.allIdxFin        = auRes.allIdxFin;
res.pertIdxFin       = auRes.pertIdxFin;
res.contIdxFin       = auRes.contIdxFin;
res.numTrialsFin     = auRes.numTrialsFin;
res.numContTrialsFin = auRes.numContTrialsFin;
res.numPertTrialsFin = auRes.numPertTrialsFin;
res.pertTrigsFin     = auRes.pertTrigsFin;

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
res.limitsAMH        = auRes.limitsAMH;        % Limits Audio Corrected for MicXHead

%Inflation Response
res.respVar      = auRes.respVar;
res.respVarM     = auRes.respVarM;
res.respVarSD    = auRes.respVarSD;
res.InflaStimVar = auRes.InflaStimVar;

%Some Final Output time series data
res.audioMinc = auRes.audioMinc;
res.audioHinc = auRes.audioHinc;
res.AuNiDelaysinc = auRes.AuNiDelaysinc;

%Check the Ni trials against the Au trials
presInd = [];
if res.numPertTrialsFin < res.numPertTrialsNi
    setPertTrials = res.allIdxFin(res.pertIdxFin);
    for ii = 1:length(res.pertIdxNi)

        ind = find(setPertTrials == res.pertIdxNi(ii));
        if ~isempty(ind)
            presInd = cat(1, presInd, ii);
        end
    end
    
    res.sensorPsv    = res.sensorP(:, presInd);
    res.sensorPSecsv = res.sensorPSec(:, presInd, :);
    res.lagTimePsv   = res.lagTimeP(presInd, :);
    res.riseTimePsv  = res.riseTimeP(presInd);
    res.OnOfValPsv   = res.OnOfValP(presInd, :);
else
    res.sensorPsv    = res.sensorP;
    res.sensorPSecsv = res.sensorPSec;
    res.lagTimePsv   = res.lagTimeP;
    res.riseTimePsv  = res.riseTimeP;
    res.OnOfValPsv   = res.OnOfValP;   
end
end

function [DRF, pF, iRF] = preAnalysisCheck(DRF, f0b)
% preAnalysisCheck(DRF) returns flags for different analyses to be 
% performed. expType will be the name of the experiment, and will likely be
% either 'Somatosensory Perturbation_Perceptual' or 'Auditory
% Perturbation_Perceptual. 
%
% This returns answers for pF (Pressure Flag) and iRF (Inflation Response
% Flag)

expType = DRF.expParam.expType;

if strcmp(expType(1:4), 'Soma') == 1
    % Laryngeal perturbations. AKA, we want to analyze the pressure, and 
    % the behavioral response to the inflation!
    pF  = 1;      %Pressure Analysis Flag
    iRF = 1;      %Inflation Response Flag
else
    % No laryngeal perturbations
    pF  = 0;      %Pressure Analysis Flag
    iRF = 0;      %Inflation Response Flag
end

if ~isfield(DRF.expParam, 'age')
    DRF.expParam.age = NaN;
end

if ~isfield(DRF.expParam, 'f0b')
    DRF.expParam.f0b = f0b;
end
end

function ss = createSortingStruct()

ss.trialType = [];
ss.sigs      = [];
ss.trigs     = [];
ss.rawData   = [];
ss.DAQin     = [];
end

function [DRFe, ss] = ammendSortingStruc(DRF, ss, j)

DRFe.dirs     = DRF.dirs;
DRFe.expParam = DRF.expParam;
DRFe.p        = DRF.p;
DRFe.audStimP = DRF.audStimP;

ss.trialType = cat(2, ss.trialType, DRF.expParam.trialType(j));
ss.sigs      = cat(2, ss.sigs, DRF.expParam.sigs(:,j));
ss.trigs     = cat(1, ss.trigs, DRF.expParam.trigs(j, :, :));
ss.rawData   = cat(1, ss.rawData, DRF.rawData);
ss.DAQin     = cat(3, ss.DAQin, DRF.DAQin);
end