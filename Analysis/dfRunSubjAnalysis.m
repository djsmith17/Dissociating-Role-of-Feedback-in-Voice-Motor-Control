function dfRunSubjAnalysis()
% dfRunSubjAnalysis() is my main script for analyzing recorded audio files 
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
AVar.participants  = {'Pilot29'};       %    List of multiple participants.
AVar.numPart       = length(AVar.participants);
AVar.runs          = {'DS1', 'DS2', 'DS3', 'DS4', 'DS5', 'DS6'}; %    List of multiple runs.
AVar.numRuns       = length(AVar.runs);
AVar.baselineFile  = 'BV1';
AVar.debug         = 0;

dirs               = dfDirs(AVar.project);

for i = 1:AVar.numPart
    participant = AVar.participants{i};
    dirs.baselineData  = fullfile(dirs.SavData, participant, AVar.baselineFile, [participant AVar.baselineFile 'DRF.mat']); % Where to find data
        
    if exist(dirs.baselineData, 'file') == 0
        fprintf('ERROR: Could not find baseline data set at %s\n', dirs.baselineData)
        return
    else
        fprintf('Loading baseline data set for %s %s\n', participant, AVar.baselineFile)
        load(dirs.baselineData) % Returns DRF
        bV = DRF;
    end
    
    for j = 1:AVar.numRuns
        run         = AVar.runs{j};
        
        dirs.SavFileDir    = fullfile(dirs.SavData, participant, run, [participant run 'DRF.mat']);                             % Where to find data
        dirs.SavResultsDir = fullfile(dirs.Results, participant, run);                                                          % Where to save results
        dirs.InflaVarDir   = fullfile(dirs.SavData, participant, 'IV1');                                                        % Where to save results

        if exist(dirs.SavResultsDir, 'dir') == 0
            mkdir(dirs.SavResultsDir)
        end
        
        if exist(dirs.InflaVarDir, 'dir') == 0
            mkdir(dirs.InflaVarDir)
        end
        
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%\n')
        if exist(dirs.SavFileDir, 'file') == 0
            fprintf('ERROR: Could not find saved data set at %s\n', dirs.SavFileDir)
            return
        else
            fprintf('Loading saved data set for %s %s\n', participant, run)
            load(dirs.SavFileDir) % Returns DRF
        end
        
        AVar.expType = DRF.expParam.expType;
        [pF, iRF] = checkDRFExpType(AVar.expType);
        aFn = 0; aFa = 1; %Audio Analysis Flag        
        
        f0b = bV.qRes.meanf0;
        [niAn, niRes] = dfAnalysisNIDAQ(dirs, DRF.expParam, DRF.DAQin, f0b, aFn, iRF, pF);
        [auAn, auRes] = dfAnalysisAudapter(dirs, DRF.expParam, DRF.rawData, f0b, aFa, iRF, niAn);

        res = combineRes(niRes, auRes);
        
        dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [participant run 'Results' res.f0Type 'DRF.mat']);
        if AVar.debug == 0
            fprintf('\nSaving Results for %s %s\n', participant, run)
            save(dirs.SavResultsFile, 'res')
        end
        
        if iRF == 1
            saveInflationResponse(dirs, res, participant, run, AVar.debug)
        end
    end
end
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
res.svIdx         = auRes.svIdx;
res.expTrigsSv    = auRes.expTrigsSv;
res.pertIdx       = auRes.pertIdx;     % The indices of the svIdx;
res.pertTrig      = auRes.pertTrig;
res.contIdx       = auRes.contIdx;     % The indices of the svIdx;
res.contTrig      = auRes.contTrig;
res.numSaveTrials = auRes.numSaveTrials;
res.numContTrials = auRes.numContTrials;
res.numPertTrials = auRes.numPertTrials;

res.balloon         = niRes.balloon;
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

function [pF, iRF] = checkDRFExpType(expType)
% This sets different variables so that the analyzes are done differently. 

if strcmp(expType, 'Somatosensory Perturbation_Perceptual') == 1
    % The perturbatron was used and I want to measure the pressure response
    % and the behavioral inflation response!
    pF  = 1;      %Pressure Analysis Flag
    iRF = 1;      %Inflation Response Flag
else
    % No perturbatron
    pF  = 0;      %Pressure Analysis Flag
    iRF = 0;      %Inflation Response Flag
end
end

function saveInflationResponse(dirs, res, participant, run, debug)

InflaVar = res.InflaStimVar;

dirs.InflaVarFile = fullfile(dirs.InflaVarDir, [participant 'IV1' 'DRF.mat']);
if debug == 0
    fprintf('Saving Inflation Stimulus Variables for %s %s\n\n', participant, run)
    save(dirs.InflaVarFile, 'InflaVar');
end
end