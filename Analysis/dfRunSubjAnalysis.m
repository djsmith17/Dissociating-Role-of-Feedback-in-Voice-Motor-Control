function dfRunSubjAnalysis()
%Analyses the microphone data from the somatosensory perturbation
%experiment. Measures the change in f0 over each trial, and each run for a
%given participant. At the end it approximates a general response to
%inflation to be used in the auditory perturbation experiment

%Require the Signal Processing Toolbox

AVar.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
AVar.participants  = {'Pilot28'}; %List of multiple participants.
AVar.numPart       = length(AVar.participants);
AVar.runs          = {'SF1'};
AVar.numRuns       = length(AVar.runs);
AVar.baselineFile  = 'BV1';
AVar.debug         = 0;

dirs               = dfDirs(AVar.project);

for i = 1:AVar.numPart
    for j = 1:AVar.numRuns
        participant = AVar.participants{i};
        run         = AVar.runs{j};
        
        dirs.baselineData  = fullfile(dirs.RecData, participant, AVar.baselineFile, [participant AVar.baselineFile 'DRF.mat']); %Where to find data
        dirs.SavFileDir    = fullfile(dirs.RecData, participant, run, [participant run 'DRF.mat']); %Where to find data
        dirs.SavResultsDir = fullfile(dirs.Results, participant, run); %Where to save full analyzed results
        
        dirs.InflaVarDir  = fullfile(dirs.RecData, participant, 'IV1');

        if exist(dirs.baselineData, 'file') == 0
            disp('ERROR')
            return
        end
        
        if exist(dirs.InflaVarDir, 'dir') == 0
            mkdir(dirs.InflaVarDir)
        end
        
        if exist(dirs.SavResultsDir, 'dir') == 0
            mkdir(dirs.SavResultsDir)
        end

        fprintf('Loading Files for %s %s\n', participant, run)
        load(dirs.baselineData) % Returns DRF
        bV = DRF;
        load(dirs.SavFileDir)   % Returns DRF
        
        AVar.expType = DRF.expParam.expType;
        [pF, iRF] = checkDRFExpType(AVar.expType);
        aFn = 0; aFa = 1; %Audio Analysis Flag        
        
        %Initialize these so I can stop worrying about it
        niAn = []; niRes = [];
        auAn = []; auRes = [];
                
        f0b = bV.qRes.meanf0;
        [niAn, niRes] = dfAnalysisNIDAQ(dirs, DRF.expParam, DRF.DAQin, f0b, aFn, iRF, pF);
        [auAn, auRes] = dfAnalysisAudapter(dirs, DRF.expParam, DRF.rawData, f0b, aFa, iRF, niAn);

        res = combineRes(niRes, auRes);
        
        dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [participant run 'Results' res.f0Type 'DRF.mat']);
        if AVar.debug == 0
            fprintf('Saving Results for %s %s\n', participant, run)
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
    fprintf('Saving Inflation Stimulus Variables for %s %s\n', participant, run)
    save(dirs.InflaVarFile, 'InflaVar');
end
end