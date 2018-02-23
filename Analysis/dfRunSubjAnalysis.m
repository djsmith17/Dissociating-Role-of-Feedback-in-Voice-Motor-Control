function dfRunSubjAnalysis()
%Analyses the microphone data from the somatosensory perturbation
%experiment. Measures the change in f0 over each trial, and each run for a
%given participant. At the end it approximates a general response to
%inflation to be used in the auditory perturbation experiment

%Require the Signal Processing Toolbox

AVar.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
AVar.participants  = {'Pilot28'}; %List of multiple participants.
AVar.numPart       = length(AVar.participants);
AVar.runs          = {'LDDDiag1'};
AVar.numRuns       = length(AVar.runs);
AVar.baselineFile  = 'BV2';
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
        load(dirs.baselineData) % Expect DRF
        bV = DRF;
        load(dirs.SavFileDir)   % Expect DRF
        
        AVar.expType = DRF.expParam.expType;
        [pF, iRF] = checkDRFExpType(AVar.expType);
        AudFlag = 1;
        
        %Initialize these so I can stop worrying about it
        niAn = []; niRes = [];
        auAn = []; auRes = [];
                
        f0b = 215; %bV.qRes.meanf0;
        [niAn, niRes] = dfAnalysisNIDAQ(dirs, DRF.expParam, DRF.DAQin, f0b, 0, pF, iRF);
        [auAn, auRes] = dfAnalysisAudapter(dirs, DRF.expParam, DRF.rawData, niAn, f0b, AudFlag);

        res = combineRes(niRes, auRes);
        
        dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [participant run 'ResultsDRF.mat']);
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

res.timef0        = auRes.timef0;
res.f0b           = auRes.f0b;

res.numContTrialsPP = auRes.numContTrialsPP;
res.numPertTrialsPP = auRes.numPertTrialsPP;
res.pertTrigPP      = auRes.pertTrigPP;

%Full Individual Trials: Mic/Head f0 Trace 
res.audioMf0TrialPert = auRes.audioMf0_pPP;
res.audioMf0TrialCont = auRes.audioMf0_cPP;
res.audioHf0TrialPert = auRes.audioHf0_pPP;
res.audioHf0TrialCont = auRes.audioHf0_cPP;
res.limitsA           = auRes.limitsA;
res.limitsAudRes      = auRes.limitsAudRes;

%Sections Trials: Mic/Head f0
res.secTime          = auRes.secTime;
res.audioMf0SecPert  = auRes.audioMf0_Secp;
res.audioMf0SecCont  = auRes.audioMf0_Secc;
res.audioHf0SecPert  = auRes.audioHf0_Secp;
res.audioHf0SecCont  = auRes.audioHf0_Secc;

%Mean Sectioned Trials: Mic/Head f0 Trace 
res.audioMf0MeanPert = auRes.audioMf0_meanp; % [MeanSigOn 90%CI MeanSigOff 90%CI]
res.audioMf0MeanCont = auRes.audioMf0_meanc;
res.audioHf0MeanPert = auRes.audioHf0_meanp;
res.audioHf0MeanCont = auRes.audioHf0_meanc;
res.limitsAmean      = auRes.limitsAmean;
res.limitsAMH        = auRes.limitsAMH;      % Limits Audio Corrected for MicXHead

%Inflation Response
res.respVar      = auRes.respVar;
res.respVarM     = auRes.respVarM;
res.respVarSD    = auRes.respVarSD;
res.InflaStimVar = auRes.InflaStimVar;
end

function [pF, iRF] = checkDRFExpType(expType)
% This sets different variables so that the analyzes are done differently. 

if strcmp(expType, 'Somatosensory Perturbation_Perceptual') == 1
    pF  = 1;      %Pressure Analysis Flag
    iRF = 1;      %Inflation Response Flag
else
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