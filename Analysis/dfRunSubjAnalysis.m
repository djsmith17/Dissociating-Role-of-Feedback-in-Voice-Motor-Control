function dfRunSubjAnalysis()
% dfRunSubjAnalysis() is my main script for analyzing recorded audio files 
% and sensor information from experiments studying voice motor control. 
% This function is set up to analyze multiple subject and runs in an 
% identical fashion
%
% One of the primary outcome measures of this analysis is the time-series
% analysis of a speaker's f0 as they complete either a somatosensory or
% auditory feedback perturbation task.
%
% This makes use of the following functions:
% -dfAnalysisNIDAQ.m
% -dfAnalysisAudapter.m
%
% See below for the following sub-functions:
% -preAnalysisCheck
% -combineRes
%
% Requires the Signal Processing Toolbox

close all
AVar.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
AVar.participants  = {'DRF1',...
                      'DRF2',...
                      'DRF4',...
                      'DRF5',...
                      'DRF6',...
                      'DRF7',...
                      'DRF8',...
                      'DRF9',...
                      'DRF10',...
                      'DRF12',...
                      'DRF13',...
                      'DRF14',...
                      'DRF15',...
                      'DRF16',...
                      'DRF17',...
                      'DRF18',...
                      'DRF19',...
                      'DRF20'}; % List of multiple participants.
AVar.numPart       = length(AVar.participants);
AVar.runs          = {'SF1', 'SF2', 'SF3', 'SF4'}; %    List of multiple runs.
AVar.numRuns       = length(AVar.runs);
AVar.baselineFile  = 'BV1';            % Baseline Voice information
AVar.debug         = 0;

dirs               = dfDirs(AVar.project);
dirs.LoadData      = dirs.SavData;

for i = 1:AVar.numPart
    participant = AVar.participants{i};
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
    for j = 1:AVar.numRuns
        run         = AVar.runs{j};
        
        % Define where to load raw data and save analyzed results
        dirs.SavFileDir    = fullfile(dirs.LoadData, participant, run, [participant run 'DRF.mat']);  % Where to find data
        dirs.SavResultsDir = fullfile(dirs.Results, participant, run);                                % Where to save results
        
        % Make sure there is a place to save results
        if exist(dirs.SavResultsDir, 'dir') == 0
            mkdir(dirs.SavResultsDir)
        end
        
        % Look for the recording sessoin raw data for this participant, then load it
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%\n')
        if exist(dirs.SavFileDir, 'file') == 0
            fprintf('ERROR: Could not find saved data set at %s\n', dirs.SavFileDir)
            return
        else
            tic
            fprintf('Loading saved data set for %s %s\n', participant, run)
            load(dirs.SavFileDir) % Returns DRF
        end
        
        % Identify the type of experiment and decide what types of analyzes
        % we need. pF: Pressure Flag; iRF: Inflation Response Flag        
        [DRF, pF, aDF] = preAnalysisCheck(DRF, f0b);
        audioFlagN = 0; % Audio Analysis Flag
        audioFlagA = 1; % Audio Analysis Flag
        f0CalcF    = 0;
        
        close all
        % Analysis on the NIDAQ raw data
        [niAn, niRes] = dfAnalysisNIDAQ(dirs, DRF.expParam, DRF.DAQin, audioFlagN, aDF, pF);
        % Analysis on the Audapter raw data
        [auAn, auRes] = dfAnalysisAudapter(dirs, DRF.expParam, DRF.rawData, niAn, audioFlagA, aDF, f0CalcF);

        % Combine Audapter and NIDAQ results into one neat MATLAB structure
        res = combineRes(niRes, auRes);
        
        dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [participant run 'ResultsDRF.mat']);
        if AVar.debug == 0
            % Save the results of this recording session
            fprintf('\nSaving Results for %s %s\n', participant, run)
            save(dirs.SavResultsFile, 'res')
            fprintf('Time Elapsed: %0.2f s\n', toc)
        end
    end
end
end

function res = combineRes(niRes, auRes)

% Identifying Subject Information
res.subject = auRes.subject;
res.run     = auRes.run;
res.curSess = auRes.curSess;
res.gender  = auRes.gender;
res.age     = auRes.age;

% Identifying Experimental Settings
res.expType  = auRes.expType;
res.AudFB    = auRes.AudFB;
res.SeqAudFB = auRes.SeqAudFB;

% Identifying Analysis Settings
res.f0Type  = auRes.f0Type;
res.etMH    = auRes.etMH;

% Raw Recorded data
res.sRate               = auRes.sRate;
res.numTrial            = auRes.numTrial;   % Total trials recorded
res.removedTrialTracker = auRes.removedTrialTracker;
res.incTrialInfo        = auRes.incTrialInfo;
res.allAuMHDelays       = auRes.allAuMHDelays; % Vector of the delays between the NIDAQ and Audapter microphone recordings
res.allAuNiDelays       = auRes.allAuNiDelays; % Vector of the delays between the NIDAQ and Audapter microphone recordings
res.prePertVoicingTimes = auRes.prePertVoicingTimes;

%NIDAQ RESULTS Collection
res.balloon         = niRes.balloon;
res.numPertTrialsNi = niRes.numPertTrials;
res.numContTrialsNi = niRes.numContTrials;
res.pertIdxNi       = niRes.pertIdx;

res.timeS           = niRes.timeS;
res.sensorP         = niRes.sensorP; % Individual Processed Perturbed trials. 
res.presSD          = niRes.presSD;  % Sensor Dynamics Structure
res.limitsP         = niRes.limitsP; % Limits for collection of individual trials

% Limits for Aligned and Meaned Pressure Recordings
res.limitsPAl   = niRes.limitsPAl;
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

% Dynamics of the Participant's Vocal Response
% This can be the response to the laryngeal pert
% This can be the response to the auditory pert
res.audioDynamics = auRes.audioDynamics;

%NIAu Delay
res.audioMinc     = auRes.audioMinc;
res.audioHinc     = auRes.audioHinc;
res.AuMHDelaysinc = auRes.AuMHDelaysinc;
res.AuNiDelaysinc = auRes.AuNiDelaysinc;
res.prePertVoicingTimeinc = auRes.prePertVoicingTimeinc;


%Result variables to be plotted after accounting for tossed trials
res.sensorPsv    = res.sensorP;
res.presSDsv     = res.presSD;

%Check the saved Au trials against the Ni trials
presInd = [];
if res.numPertTrialsFin < res.numPertTrialsNi && strcmp(res.expType, 'Somatosensory Perturbation_Perceptual')
    setPertTrials = res.allIdxFin(res.pertIdxFin);
    for ii = 1:length(res.pertIdxNi)

        ind = find(setPertTrials == res.pertIdxNi(ii));
        if ~isempty(ind)
            presInd = cat(1, presInd, ii);
        end
    end
    
    res.sensorPsv          = res.sensorPsv(:, presInd);
    res.presSDsv.sensorSec = res.presSDsv.sensorSec(:, presInd, :);
    res.presSDsv.lagTimes  = res.presSDsv.lagTimes(presInd, :);
    res.presSDsv.riseTimes = res.presSDsv.riseTimes(presInd);
    res.presSDsv.OnOffVal  = res.presSDsv.OnOffVal(presInd, :);
end
end

function [DRF, pF, aDF] = preAnalysisCheck(DRF, f0b)
% preAnalysisCheck(DRF) returns flags for different analyses to be 
% performed. expType will be the name of the experiment, and will likely be
% either 'Somatosensory Perturbation_Perceptual' or 'Auditory
% Perturbation_Perceptual. 
%
% This returns
% -pF (Pressure Flag)
% -aDF (Audio Dynamics Flag)

expType = DRF.expParam.expType;
if strcmp(expType, 'Somatosensory Perturbation_Perceptual') == 1
    % Laryngeal perturbations. AKA, we want to analyze the pressure, and 
    % the behavioral response to the inflation!
    pF  = 1;      %Pressure Analysis Flag
    aDF = 1;      %Audio Dynamics Flag
elseif strcmp(expType, 'Auditory Perturbation_Perceptual') == 1
    % Auditory perturbations. AKA, no pressure, but also behavioral response to the pitch shifted stimuli
    pF  = 0;      %Pressure Analysis Flag
    aDF = 2;      %Audio Dynamics Flag
else
    % No laryngeal perturbations
    pF  = 0;      %Pressure Analysis Flag
    aDF = 0;      %Audio Dynamics Flag
end

if ~isfield(DRF.expParam, 'age')
    DRF.expParam.age = NaN;
end

if ~isfield(DRF.expParam, 'f0b')
    DRF.expParam.f0b = f0b;
end
end