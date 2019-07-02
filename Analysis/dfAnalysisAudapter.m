function [auAn, auRes] = dfAnalysisAudapter(dirs, expParam, rawData, niAn, AudFlag, aDF, f0CalcF)
% [auAn, auRes] = dfAnalysisAudapter(dirs, expParam, rawData, f0b, AudFlag, iRF, niAn)
% This function analyzes the raw audio data that was recorded by Audapter 
% in the experiments measuring changes in f0. It first does a
% preprocessing step where it identifies any temporal errors in
% production, and also identifies and corrects for any lags in recording.
% Once all the data are processed, they are handed off
% to a function to do the frequency analysis of the audio signals. 
%
% An important distinction in the processing step is the difference between
% the recorded trials that are to be saved and further analyzed, and those
% that will be thrown out. Any variables that refer to the full trial set
% will not have any additional prefix/suffix. Variables refering to the 
% set of trials that passed the 'temporal' preprocessing and will be 
% saved for further analyses will have a prefix/suffix of 'Svt'.
% 
% dirs:     The set of directories we are currently working in 
% expParam: The experimental parameters of the recorded experiment
% rawData:  Raw Audapter data structures
% f0b:      Baseline fundamental frequency, recorded from baseline trials
% AudFlag:  Flag to check if analyses of audio data should be performed
% iRF:      Inflation Response Flag; should the inflation response be calculated
% niAn:     Analysis variables used to analyze NIDAQ data
%
% auAn:  Analysis variables used to analyze Audapter data
% auRes: Structure of result vars that are needed for stats and plotting
%
% This function calls the following functions
% dfAnalysisAudio.m
%
% Requires the Signal Processing Toolbox

% Identify some Experimental variables
auAn = setAudapterAnalysisParams();

auAn.expType   = expParam.expType;
auAn.subject   = expParam.subject;
auAn.run       = expParam.run;
auAn.curSess   = expParam.curSess;
auAn.gender    = expParam.gender;
auAn.f0AnaFile = [auAn.subject auAn.run 'f0Analysis.mat'];
auAn.age       = expParam.age;
auAn.f0b       = expParam.f0b;
auAn.AudFB     = expParam.AudFB;
auAn.AudFBSw   = expParam.AudFBSw;

fprintf('\nStarting Analysis for %s, %s with f0 of %0.2f Hz\n', auAn.subject, auAn.run, auAn.f0b)

% Idenitfy some Recording Variables
auAn.sRate    = expParam.sRateAnal;
auAn.sRateNi  = niAn.sRate;
auAn.frameLen = expParam.frameLenDown;
auAn.trialLen = expParam.trialLen;
auAn.numSamp  = auAn.sRate*auAn.trialLen;

auAn.numTrial   = expParam.numTrial;
auAn.trialType  = expParam.trialType;
auAn.expTrigs   = expParam.trigs(:,:,1);
auAn.anaTrigs   = expParam.trigs(:,:,3);

if isfield(expParam, 'incTrialInfo')
    auAn.incTrialInfo = expParam.incTrialInfo;
end

% Is there a field for SEQAudFB? Sometimes AudFB changes between trials
if isfield(expParam, 'SeqAudFB')
    auAn.SeqAudFB   = expParam.SeqAudFB;
    auAn.SeqAudFBSw = expParam.SeqAudFBSw;
else
    auAn.SeqAudFB      = cell(1, auAn.numTrial);
    [auAn.SeqAudFB{:}] = deal(auAn.AudFB);
    auAn.SeqAudFBSw    = repmat(auAn.AudFBSw, 1,  auAn.numTrial);
end

for ii = 1:auAn.numTrial
    data = rawData(ii);       % Get the data from this trial
    
    trialVar.rawMic  = data.signalIn;                    % Microphone
    trialVar.rawHead = data.signalOut;                   % Headphones
    trialVar.rms     = data.rms(:,1);                    % RMS recording
    trialVar.auTrigs = auAn.anaTrigs(ii,:);              % Perturbation Triggers (in Audapter points)
    trialVar.AudFB   = auAn.SeqAudFBSw(ii);              % Auditory Feedback Used
    trialVar.typeIdx = auAn.trialType(ii);               % Trial Type 0(Control), 1 (Perturbed)
    trialVar.type    = auAn.types{trialVar.typeIdx + 1}; % Trial Type (Words)
    
    trialVar = setNIDAQTrialVars(niAn, auAn, ii, trialVar);
    
    % Preprocessing step identifies time-series errors in recording/vocalization
    % Making use of MH class object to handle the signal processing. See
    % that script for more information
    MH = MicHeadAlignProcess(auAn, trialVar);
    
    % Save all trials (and delay calcs), regardless of eventual exclusion
    auAn.audioM              = cat(2, auAn.audioM, MH.processedMic);
    auAn.audioH              = cat(2, auAn.audioH, MH.processedHead);
    auAn.prePertVoicingTimes = cat(1, auAn.prePertVoicingTimes, MH.prePertVoicingTime);
    auAn.allAuMHDelays       = cat(1, auAn.allAuMHDelays, MH.AuMHDelay);
    auAn.allAuNiDelays       = cat(1, auAn.allAuNiDelays, MH.AuNIDelay);
    
    % Identify if trial should be tossed
    if MH.saveT == 0     % Don't save the trial :(
        fprintf('%s Trial %d (%s) excluded due to %s\n', auAn.curSess, ii, trialVar.type, MH.saveTmsg)
        removedTrial = {['Trial ' num2str(ii)], MH.saveTmsg};
        auAn.removedTrialTracker = cat(1, auAn.removedTrialTracker, removedTrial);
        
    elseif MH.saveT == 1 % Save the Trial        
        auAn.allIdxPreProc = cat(1, auAn.allIdxPreProc, ii); % Save the experimental index
    end
end

% The Audio Analysis
auAn = dfAnalysisAudio(dirs, auAn, AudFlag, aDF, f0CalcF);

% (inc)luded trials following Audio Analysis trial removal
auAn.audioMinc             = auAn.audioM(:, auAn.svf0Idx);
auAn.audioHinc             = auAn.audioH(:, auAn.svf0Idx);
auAn.AuMHDelaysinc         = auAn.allAuMHDelays(auAn.svf0Idx);
auAn.AuNiDelaysinc         = auAn.allAuNiDelays(auAn.svf0Idx);
auAn.prePertVoicingTimeinc = auAn.prePertVoicingTimes(auAn.svf0Idx);

lims  = dfIdentifyAudioLims(auAn);
auRes = packResults(auAn, lims);
end

function auAn = setAudapterAnalysisParams()

auAn.AnaType   = 'Audapter';% Analysis of data recorded with Audapter
auAn.expType   = [];        % Somatosensory vs Auditory
auAn.subject   = [];        % Participant name
auAn.run       = [];        % Run name
auAn.curSess   = [];        % Combination of Participant name and Run name
auAn.f0Type    = 'Praat';   % Which algorithm for-pitch tracking
auAn.f0AnaFile = [];        % Where we will save the f0 trace analyzed data
auAn.gender    = [];        % Participant gender
auAn.age       = [];        % Participant age
auAn.f0b       = [];        % Participant f0 from baseline recording
auAn.AudFB     = [];        % Form of Auditory Feedback used in experimental run
auAn.AudFBSw   = [];        % Auditory Feedback settings used in experimental run
auAn.SeqAudFB  = [];        % If different AudFB was presented from trial-to-trial, track it here. 
auAn.SeqAudFBSw = [];       % If different AudFB was presented from trial-to-trial, track it here. 

auAn.sRate     = []; % Sampling Rate of Audapter (down-sampled)
auAn.sRateNi   = []; % Sampling Rate of NIDAQ
auAn.frameLen  = []; % Frame Length of Audapter recording (down-sampled)
auAn.trialLen  = []; % Length of recording (s)
auAn.numSamp   = []; % Length of recording (points)
auAn.rmsThresh = []; % rms level for voicing

auAn.time      = []; % Time Vector based on numSamp
auAn.audioM    = []; % Mic data that has received temporal preprocessing (all recorded trials)
auAn.audioH    = []; % Head data that has received temporal preprocessing (all recorded trials)
auAn.numTrial  = []; % Number of trials recorded
auAn.trialType = []; % Key for identifying Control (0) & Perturbed (1) trials
auAn.types     = {'Control', 'Perturbed'};
auAn.expTrigs  = []; % Trigger Onset and Offset (Time) (all recorded trials)
auAn.anaTrigs  = []; % Trigger Onset and Offset (Points; Audadapter) 
auAn.removedTrialTracker = {}; % List of Trials that were automatically thrown out during Analysis
auAn.incTrialInfo = []; % Record of manual trial removal
auAn.allAuMHDelays = []; % Vector of the delays between the Audapter microphone and Headphone recordings
auAn.allAuNiDelays = []; % Vector of the delays between the NIDAQ and Audapter microphone recordings
auAn.prePertVoicingTimes = []; % Vector of the times that are recorded of the participant vocalizing pre-perturbation

auAn.allIdxPreProc = []; % Vector of indicies of all recorded trials saved for further analyses.
auAn.audioMSvt     = []; % Microphone recordings for the trials saved for further analyses
auAn.audioHSvt     = []; % Headphone recordings for the trials saved for further analyses
auAn.numTrialSvt   = []; % Number of trials saved for further analyses
auAn.trialTypeSvt  = []; % Key for identifying Control (0) & Perturbed (1) trials
auAn.expTrigsSvt   = []; % Trigger Onset and Offset (Time) for trials saved for further analyses
end

function trialVar = setNIDAQTrialVars(niAn, auAn, ii, trialVar)
% trialVar = setNIDAQTrialVars(niAn, auAn, ii, trialVar) is a simple
% function to set up the fields needed to describe the nidaq behavior.
% Sometimes NIDAQ data was not recorded, and we want the MHAlignProcess to
% still function correctly. This sets default values in that case.

if isfield(niAn, 'audioM')
    trialVar.fsNI        = niAn.sRate;
    trialVar.trialTimeNI = niAn.trialLen;
    trialVar.numSampNI   = niAn.numSamp;
    trialVar.timeNI      = niAn.time;
    trialVar.rawMicNI    = niAn.audioM(:,ii);          % Microphone (NIDAQ)
    trialVar.pressureNI  = niAn.sensorPz(:,ii);
    trialVar.expTrigsNI  = niAn.expTrigs(ii, :);
    
    if ismember(ii, niAn.pertIdx) && ~isempty(niAn.presSD.TrigTime)
        thisIdx = find(ii == niAn.pertIdx);
        trialVar.pressureTrigs = niAn.presSD.TrigTime(thisIdx, :);
        trialVar.presLagTimes  = niAn.presSD.lagTimes(thisIdx, :);
        trialVar.presRiseTimes = niAn.presSD.riseTimes(thisIdx, :);
    else
        trialVar.pressureTrigs = trialVar.expTrigsNI;
        trialVar.presLagTimes  = [0 0];
        trialVar.presRiseTimes = [0 0];
    end
else
    trialVar.fsNI        = 8000;
    trialVar.trialTimeNI = 4;
    trialVar.numSampNI   = trialVar.trialTimeNI*trialVar.fsNI;
    trialVar.timeNI      = linspace(0, trialVar.trialTimeNI, trialVar.numSampNI);
    trialVar.rawMicNI    = resample(trialVar.rawMic, trialVar.fsNI, auAn.sRate);
    trialVar.pressureNI  = zeros(size(trialVar.timeNI));
    trialVar.expTrigsNI    = [0 0];
    trialVar.pressureTrigs = [0 0];
    trialVar.presLagTimes  = [0 0];
    trialVar.presRiseTimes = [0 0];
end
end

function res = packResults(auAn, lims)
% res = packResults(auAn, lims) takes the values stored in the analysis 
% variables structure and repackages the important variables into a new 
% structure that have common names between other analysis methods. 
% This makes it easy to switch back and forth between different result 
% structures for plotting, and it makes the plotting functions reliant on
% a a uniform naming style, that does not change, even when the analysis
% methods/names may.

% Information about the experiment/subject
res.subject = auAn.subject;
res.run     = auAn.run;
res.curSess = auAn.curSess;
res.gender  = auAn.gender;
res.age     = auAn.age;

res.expType  = auAn.expType;          % Somatosensory Perturbation_Perceptual, etc
res.AudFB    = auAn.AudFB;            % Voice Feedback, Masking Noise
res.SeqAudFB = auAn.SeqAudFB;         % Same as above, but keeps track of individual trial differences

res.f0Type = auAn.f0Type;             % Type of pitch-shift perturbation used
res.etMH   = auAn.etMH;               % Time it took to run the f0 analysis

% Raw Recorded data
res.sRate         = auAn.sRate;
res.numTrial      = auAn.numTrial;    % Total trials recorded
res.audioM        = auAn.audioM;      % Raw microphone recording (no lag correction)
res.audioH        = auAn.audioH;      % Raw headphone recording  (no lag correction)
res.allAuMHDelays = auAn.allAuMHDelays; % Vector of the delays between the NIDAQ and Audapter microphone recordings
res.allAuNiDelays = auAn.allAuNiDelays; % Vector of the delays between the NIDAQ and Audapter microphone recordings
res.prePertVoicingTimes = auAn.prePertVoicingTimes;

res.removedTrialTracker = auAn.removedTrialTracker; % Result of Automatic Trial Exclusion
res.incTrialInfo = auAn.incTrialInfo; % Result of Manual Trial Exclusion

% Post temporal processing data
res.allIdxPreProc = auAn.allIdxPreProc; % Vector of indicies of recorded trials saved (post temporal processing)

% Audio f0 analysis
res.timef0        = auAn.timef0;
res.f0b           = auAn.trialf0M;

res.svf0Idx      = auAn.svf0Idx;
res.expTrigsf0Sv = auAn.expTrigsf0Sv;
res.pertf0Idx    = auAn.pertf0Idx;
res.contf0Idx    = auAn.contf0Idx;

res.numTrialsPP     = auAn.numTrialsPP;
res.numContTrialsPP = auAn.numContTrialsPP;
res.numPertTrialsPP = auAn.numPertTrialsPP;
res.pertTrigPP      = auAn.pertTrigsR;

% Which trials did I end up saving and plotting at the end of the day?
res.allIdxFin        = res.svf0Idx;
res.pertIdxFin       = res.pertf0Idx;
res.contIdxFin       = res.contf0Idx;
res.numTrialsFin     = res.numTrialsPP;
res.numContTrialsFin = res.numContTrialsPP;
res.numPertTrialsFin = res.numPertTrialsPP;
res.pertTrigsFin     = res.pertTrigPP;

%Individual Full Trials: Mic/Head f0 Trace 
res.audioMf0TrialPert = auAn.audioMf0p;
res.audioMf0TrialCont = auAn.audioMf0c;
res.audioHf0TrialPert = auAn.audioHf0p;
res.audioHf0TrialCont = auAn.audioHf0c;
res.limitsA           = lims.audio;

%Individual Sectioned Trials: Mic/Head f0 Trace
res.secTime          = auAn.secTime;
res.audioMf0SecPert  = auAn.audioMf0_Secp;
res.audioMf0SecCont  = auAn.audioMf0_Secc;
res.audioHf0SecPert  = auAn.audioHf0_Secp;
res.audioHf0SecCont  = auAn.audioHf0_Secc;

%Mean Sectioned Trials: Mic/Head f0 Trace 
res.audioMf0MeanPert = auAn.audioMf0_meanp; % [MeanSigOn 90%CI MeanSigOff 90%CI]
res.audioMf0MeanCont = auAn.audioMf0_meanc;
res.audioHf0MeanPert = auAn.audioHf0_meanp;
res.audioHf0MeanCont = auAn.audioHf0_meanc;
res.limitsAmean      = lims.audioMean;

% Dynamics of the Participant's Vocal Response
res.audioDynamics = auAn.audioDynamics;

% Some Final Output time series data
res.audioMinc     = auAn.audioMinc;
res.audioHinc     = auAn.audioHinc;
res.AuMHDelaysinc = auAn.AuMHDelaysinc;
res.AuNiDelaysinc = auAn.AuNiDelaysinc;
res.prePertVoicingTimeinc = auAn.prePertVoicingTimeinc;
end