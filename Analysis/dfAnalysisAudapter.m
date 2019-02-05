function [auAn, auRes] = dfAnalysisAudapter(dirs, expParam, rawData, f0b, AudFlag, iRF, niAn)
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

auAn.time       = (0:1/auAn.sRate:(auAn.numSamp-1)/auAn.sRate)';
auAn.numTrial   = expParam.numTrial;
auAn.trialType  = expParam.trialType;
auAn.expTrigs   = expParam.trigs(:,:,1);
auAn.anaTrigs   = expParam.trigs(:,:,3);

if isfield(expParam, 'incTrialInfo')
    auAn.incTrialInfo = expParam.incTrialInfo;
end

if isfield(expParam, 'SeqAudFB')
    auAn.SeqAudFB   = expParam.SeqAudFB;
    auAn.SeqAudFBSw = expParam.SeqAudFBSw;
else
    auAn.SeqAudFB      = cell(1, auAn.numTrial);
    [auAn.SeqAudFB{:}] = deal(auAn.AudFB);
    auAn.SeqAudFBSw    = repmat(auAn.AudFBSw, 1, auAn.numTrial);
end

for ii = 1:auAn.numTrial
    data = rawData(ii);       % Get the data from this trial
    
    Mraw     = data.signalIn;     % Microphone
    Hraw     = data.signalOut;    % Headphones
    rms      = data.rms(:,1);     % RMS recording
    anaTrigs = auAn.anaTrigs(ii,:);
    AudFBSw  = auAn.SeqAudFBSw(ii);
    typeIdx  = auAn.trialType(ii);
    type     = auAn.types{typeIdx + 1};
    
    if isfield(niAn, 'audioM')
        MrawNi = niAn.audioM(:,ii);          % Microphone (NIDAQ) 
    else
        MrawNi = resample(Mraw, 8000, auAn.sRate);
    end
    
    % Preprocessing step identifies time-series errors in production/recording
    [mic, head, preProSt] = preProcAudio(auAn, Mraw, Hraw, rms, MrawNi, anaTrigs, AudFBSw);
    
    auAn.audioM = cat(2, auAn.audioM, mic);  % Save all trials, regardless of eventual exclusion
    auAn.audioH = cat(2, auAn.audioH, head); % Save all trials, regardless of eventual exclusion

    auAn.allAuMHDelays = cat(1, auAn.allAuMHDelays, preProSt.AuMHdelay);
    auAn.allAuNiDelays = cat(1, auAn.allAuNiDelays, preProSt.AuNidelay);
    if preProSt.saveT == 0     % Don't save the trial :(
        fprintf('%s Trial %d (%s) excluded due to %s\n', auAn.curSess, ii, type, preProSt.saveTmsg)
        removedTrial = {['Trial ' num2str(ii)], preProSt.saveTmsg};
        auAn.removedTrialTracker = cat(1, auAn.removedTrialTracker, removedTrial);
        
    elseif preProSt.saveT == 1 % Save the Trial        
        auAn.allIdxPreProc = cat(1, auAn.allIdxPreProc, ii); % Save the experimental index
    end
end

% The Audio Analysis
f0Flag = 1;
auAn = dfAnalysisAudio(dirs, auAn, AudFlag, iRF, f0Flag);

auAn.audioMinc = auAn.audioM(:, auAn.svf0Idx);
auAn.audioHinc = auAn.audioH(:, auAn.svf0Idx);
auAn.AuMHDelaysinc = auAn.allAuMHDelays(auAn.svf0Idx);
auAn.AuNiDelaysinc = auAn.allAuNiDelays(auAn.svf0Idx);

lims  = identifyLimits(auAn);
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

auAn.allIdxPreProc = []; % Vector of indicies of all recorded trials saved for further analyses.
auAn.audioMSvt     = []; % Microphone recordings for the trials saved for further analyses
auAn.audioHSvt     = []; % Headphone recordings for the trials saved for further analyses
auAn.numTrialSvt   = []; % Number of trials saved for further analyses
auAn.trialTypeSvt  = []; % Key for identifying Control (0) & Perturbed (1) trials
auAn.expTrigsSvt   = []; % Trigger Onset and Offset (Time) for trials saved for further analyses
end

function [micP, headP, pp] = preProcAudio(An, micR, headR, rms, micRNi, auTrigs, AudFBSw)
% [micP, headP, AuNidelay, pp] = preProcAudio(An, micR, headR, micRNi, auTrigs)
% This function performs preprocessing on the time-series recorded audio 
% data before frequency analysis methods are applied. This identifies
% delays between the Audapter recorded audio and the NIDAQ recorded audio,
% and processing delays between the Audapter microphone and headphone
% recordings. The Audapter audio signals are shifted after the delays are
% indentified.
%
% This script also calculates the time of voice onset and identifies if 
% the participant started too late (during the pre-perturbation period), 
% or had a voice break If either of these are the case, the trial is thrown
% out for further analyses.
% 
% Inputs:
% An:      Analysis variables structure
% micR:    Raw Microphone signal (Audapter)
% headR:   Raw Headphone signal  (Audapter)
% micRNI:  Raw Microphone signal (NIDAQ)
% auTrigs: Trigger points (Au) of perturbation onset and offset (per trial) 
%
% Outputs:
% micP:      Processed Microphone signal
% headP:     Processed Headphone signal
% pp:        Preprocessing results structure. This has information
%            regarding the envelope of the recorded audio file, and 
%            if the participant started late, or has a voice break. 

pp.expType     = An.expType;
pp.rawMic      = double(micR);  % Convert to data type double
pp.rawHead     = double(headR); % Convert to data type double
pp.rms         = double(rms);   % Convert to data type double
pp.fs          = An.sRate;      % Sampling rate (Audapter)
pp.frameLen    = An.frameLen;   % Frame rate of recording (After downsampling)
pp.trialLen    = length(pp.rawMic);
pp.trialTime   = pp.trialLen/pp.fs;
pp.t           = linspace(0, pp.trialTime, pp.trialLen);
pp.auTrigs     = auTrigs;
pp.AudFBSw     = AudFBSw;
pp.frameDel    = 7;
pp.rmsThresh   = 0.011;
pp.voiceOnM    = 2;

pp.micRNi      = micRNi;
pp.fsNI        = An.sRateNi;  % Sampling rate (NIDAQ)
pp.trialLenNi  = length(micRNi);
pp.trialTimeNi = pp.trialLenNi/pp.fsNI;
pp.tNi         = linspace(0, pp.trialTimeNi, pp.trialLenNi);

pp.numSamp     = pp.trialTimeNi*pp.fs;

pp.envCutOff = 40;   % Cutoff frequency for enveloping the audio
pp.thresh    = 0.30; % Threshold of Decimal amount of full peak height
pp.breakTol  = 0.1;  % Voice Break Tolerance; Time (s)

% 4th order low-pass butter filter settings
[B, A] = butter(4, pp.envCutOff/(pp.fs/2));

% Envelope the signal by low-pass filtering (change in amplitude/time ~RMS)
pp.env = filter(B, A, abs(pp.rawMic));  

% Largest peak in the envelope theoretically occurs during voicing
pp.maxPeak = max(pp.env);

% Find values that are within threshold of max 'voicing' value
pp.threshIdx = find(pp.env > pp.thresh*pp.maxPeak); 

% First index of the theoretical useable signal (Voice onset)
if pp.voiceOnM == 1
    pp.voiceOnsetInd = pp.threshIdx(1);
    pp.voiceOnsetT   = pp.t(pp.voiceOnsetInd);
else
    voicingInd = find(pp.rms > pp.rmsThresh);
    pp.voiceOnsetInd = (voicingInd(1) - pp.frameDel)*pp.frameLen;
    pp.voiceOnsetT   = pp.t(pp.voiceOnsetInd);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pp.preVOnsetTime   = 0.05; %50ms before voice onset
pp.VOnsetFrame     = floor(pp.voiceOnsetInd/pp.frameLen);
pp.preVOnsetFrames = floor(pp.preVOnsetTime*pp.fs/pp.frameLen);

pp.preVoiceRange = (-pp.preVOnsetFrames:0)+pp.VOnsetFrame;
if sum(pp.preVoiceRange <= 0) > 0
    pp.preVOnsetRMS = 0;
else
    pp.preVOnsetRMS = mean(pp.rms(pp.preVoiceRange));
end

% Find the delay between NIDAQ recording and Audapter recording
pp.micRds     = resample(pp.rawMic, pp.fsNI, pp.fs);         % Downsample the Audapter recording
pp.AuNidelay  = xCorrTimeLag(pp.micRNi, pp.micRds, pp.fsNI); % Perform xCorr between NIDAQ and Audapter. Expect that NIDAQ leads Audapter
pp.AuNidelayP = pp.AuNidelay*pp.fs;                          % Convert to points

% Adjust Triggers against NIDAQ only if we are using Laryngeal Pert Exp.
% Otherwise adjsut based on VoiceOnset, which is what Audapter does in PSR
if strcmp(pp.expType(1:3), 'Som')
    pp.adjustedDelay = pp.AuNidelayP;
    auTrigsAuNi = pp.auTrigs + pp.adjustedDelay;
else
    pp.adjustedDelay = pp.voiceOnsetInd;
    auTrigsAuNi = pp.auTrigs + pp.adjustedDelay;
end

% Aim to section audio at 0.5s pre-onset to 1.0s post-offset.
preOn   = 0.5*pp.fs;
postOff = 1.0*pp.fs;

% Audio points on either side of the perturbation period.
pp.analysisSec(1) = auTrigsAuNi(1) - preOn;   % Where to start the Analysis period
pp.analysisSec(2) = auTrigsAuNi(2) + postOff; % Where to end the Analysis period
pp.analysisPoints = pp.analysisSec(1):pp.analysisSec(2);
pp.analysisFrames = round(pp.analysisSec(1)/pp.frameLen):round(pp.analysisSec(2)/pp.frameLen);

% Check the voice onset time against when we want to start analyzing data
pp.voiceOnsetLate = pp.analysisSec(1) < pp.voiceOnsetInd;

% Check the rest of the signal following the first analysis index...are there any dead zones??
pp.fallOffLog = pp.rms(pp.analysisFrames) < pp.preVOnsetRMS;
pp.chk4Break  = sum(pp.fallOffLog) > pp.breakTol*pp.fs/pp.frameLen; % Last longer than break tolerance

% Find the delay between Audapter Headphone and Microphone
if pp.AudFBSw == 2 % No Headphone Out
    pp.AuMHdelay = (pp.frameLen*(pp.frameDel-1))/pp.fs;
else
    pp.AuMHdelay = xCorrTimeLag(pp.rawMic(pp.analysisPoints), pp.rawHead(pp.analysisPoints), pp.fs);   % Expect Mic leads Head
end
pp.AuMHdelayP = pp.AuMHdelay*pp.fs; % Convert to points

% Align the Microphone and Headphones
if pp.AuMHdelayP >= 0
    micAuAl  = pp.rawMic(1:(end-pp.AuMHdelayP));
    headAuAl = pp.rawHead((pp.AuMHdelayP+1):end);
else
    fprintf('Odd xCorr Here\n')
    micAuAl  = pp.rawMic;
    headAuAl = pp.rawHead;
end

% Adjust for delay between Audapter and NIDAQ
if pp.adjustedDelay > 0 % As long as the delay is non 0
    micAuNi  = micAuAl(pp.adjustedDelay:end);
    headAuNi = headAuAl(pp.adjustedDelay:end);
else
    micAuNi  = micAuAl;
    headAuNi = headAuAl;
end

if pp.voiceOnsetLate
    saveT    = 0;  
    saveTmsg = 'Participant started too late!!';
elseif pp.chk4Break
    saveT    = 0;
    saveTmsg = 'Participant had a voice break!!';
elseif length(micAuNi) < pp.numSamp
    saveT    = 0;
    saveTmsg = 'Recording not long enough';
    
    diffLen = pp.numSamp - length(micAuNi);
    micAuNi  = [micAuNi; zeros(diffLen, 1)];
    headAuNi = [headAuNi; zeros(diffLen, 1)];
else
    saveT    = 1;
    saveTmsg = 'Everything is good';
end

% Grab the full numSamp so they can be concatenated cleanly
micP    = micAuNi(1:pp.numSamp);
headP   = headAuNi(1:pp.numSamp);

pp.saveT    = saveT;    % Save trial or no?
pp.saveTmsg = saveTmsg; % Reason, if any the trial was thrown out
end

function [timeSet, delaySet] = MHdelayChunked(sig1, sig2, fs)

numSamp = length(sig1);
chunkL = 0.05;
chunkP = fs*chunkL;
numChunk = floor(numSamp/chunkP);

timeSet  = zeros(numChunk, 1);
delaySet = zeros(numChunk, 1);
for ii = 1:numChunk
    set = (1:chunkP) + (ii-1)*chunkP;
    
    timeChunk = set(1)/fs;
    sig1Chunk = sig1(set);
    sig2Chunk = sig2(set);
    
    delay = xCorrTimeLag(sig1Chunk, sig2Chunk, fs);
    timeSet(ii)  = timeChunk;
    delaySet(ii) = delay*1000;
end
end

function drawPreProcessDiagnostic(pp)

plotPos = [1800 10];
plotDim = [1200 900];

ppDiag = figure('Color', [1 1 1]);
set(ppDiag, 'Position', [plotPos plotDim],'PaperPositionMode','auto')

ha = tight_subplot(2,1,[0.1 0.05],[0.12 0.15],[0.08 0.08]);

axes(ha(1))
plot(pp.t, pp.rawMic)
hold on
plot(pp.t, pp.env, 'y')
hold on
plot([pp.voiceOnsetT pp.voiceOnsetT ], [-0.2 0.2])
hold on
plot([pp.t(pp.auTrigs(1)) pp.t(pp.auTrigs(1))], [-0.2 0.2], 'k--')
hold on
plot([pp.t(pp.auTrigs(2)) pp.t(pp.auTrigs(2))], [-0.2 0.2], 'k--')
box off
axis([0 6 -0.25 0.25])

axes(ha(2))
plot(pp.tNi, pp.micRNi)
title(num2str(pp.AuNidelay))
axis([0 4 -0.25 0.25])
box off

end

function timeLag = xCorrTimeLag(sig1, sig2, fs)
% xCorrTimeLag(sig1, sig2, fs) calculates the lag between two (seemingly) 
% identical time based signals. 
%
% if timeLag is positive, then sig1 leads sig2. 
% if timeLag is negative, then sig1 lags sig2.

% Simple crosscorrelation between two signals
% Finds the largest peak of the result
[r, lags]    = xcorr(sig1, sig2);
[~, peakInd] = max(r);
maxLag       = lags(peakInd);
timeLag      = maxLag/fs;
timeLag      = -timeLag;
end

function lims = identifyLimits(An)
% lims = identifyLimits(An) calculates limits of analyzed data 
% so that the limits used in plotting are dynamic and fit the data. 
%
% lims is a structure of the resultant limits [X1 X2 Y1 Y2] for given data
% This sub function is redundant within other functions and might
% eventually become its own function

%%%%%%%%%%%lims.audio%%%%%%%%%%%
%Individual Full Trials (Perturbed): f0 Audio
if ~isempty(An.audioMf0p)
    pertTrialsM = An.audioMf0p;
    pertTrialsH = An.audioHf0p;
    sec = 100:700;

    uLMa = max(pertTrialsM(sec,:));
    uLMa(find(uLMa > 200)) = 0;
    lLMa = min(pertTrialsM(sec,:));
    lLMa(find(lLMa < -300)) = 0;

    uLM = round(max(uLMa)) + 20;
    lLM = round(min(lLMa)) - 20;
       
    uLHa = max(pertTrialsH(sec,:));
    uLHa(find(uLHa > 200)) = 0;
    lLHa = min(pertTrialsH(sec,:));
    lLHa(find(lLHa < -300)) = 0;
    
    uLH = round(max(uLHa)) + 20;
    lLH = round(min(lLHa)) - 20;
    
    if uLH > uLM
        uLMH = uLH;
    else
        uLMH = uLM;
    end

    if lLH < lLM
        lLMH = lLH;
    else
        lLMH = lLM;
    end
    
    lims.audioM         = [0 4 lLM uLM];
    lims.audioAudRespMH = [0 4 lLMH uLMH];
else
    lims.audioM         = [0 4 -20 20];
    lims.audioAudRespMH = [0 4 -100 100];
end

%%%%%%%%%%%lims.audioMean%%%%%%%%%%%
%Mean Sectioned Trials (Perturbed): f0 Audio 
if ~isempty(An.audioMf0_meanp)
    [~, Imax] = max(An.audioMf0_meanp(:,1)); %Max Pert Onset
    upBoundOn = round(An.audioMf0_meanp(Imax,1) + An.audioMf0_meanp(Imax,2) + 10);
    [~, Imin] = min(An.audioMf0_meanp(:,1)); %Min Pert Onset
    lwBoundOn = round(An.audioMf0_meanp(Imin,1) - An.audioMf0_meanp(Imin,2) - 10);

    [~, Imax] = max(An.audioMf0_meanp(:,3)); %Max Pert Offset
    upBoundOf = round(An.audioMf0_meanp(Imax,3) + An.audioMf0_meanp(Imax,4) + 10);
    [~, Imin] = min(An.audioMf0_meanp(:,3)); %Min Pert Offset
    lwBoundOf = round(An.audioMf0_meanp(Imin,3) - An.audioMf0_meanp(Imin,4) - 10);

    if upBoundOn > upBoundOf
        upBoundM = upBoundOn;
    else
        upBoundM = upBoundOf;
    end

    if lwBoundOn < lwBoundOf
        lwBoundM = lwBoundOn;
    else
        lwBoundM = lwBoundOf;
    end

    lims.audioMean = [-0.5 1.0 lwBoundM upBoundM];
else
    lims.audioMean = [-0.5 1.0 -50 50];
end

%%%%%%%%%%%lims.audioMH%%%%%%%%%%%%
%Mean Sectioned Trials (Perturbed): f0 Audio 
if ~isempty(An.audioHf0_meanp)
    [~, Imax] = max(An.audioHf0_meanp(:,1)); %Max Pert Onset
    upBoundOn = round(An.audioHf0_meanp(Imax,1) + An.audioHf0_meanp(Imax,2) + 10);
    [~, Imin] = min(An.audioHf0_meanp(:,1)); %Min Pert Onset
    lwBoundOn = round(An.audioHf0_meanp(Imin,1) - An.audioHf0_meanp(Imin,2) - 10);

    [~, Imax] = max(An.audioHf0_meanp(:,3)); %Max Pert Offset
    upBoundOf = round(An.audioHf0_meanp(Imax,3) + An.audioHf0_meanp(Imax,4) + 10);
    [~, Imin] = min(An.audioHf0_meanp(:,3)); %Min Pert Offset
    lwBoundOf = round(An.audioHf0_meanp(Imin,3) - An.audioHf0_meanp(Imin,4) - 10);

    if upBoundOn > upBoundOf
        upBoundH = upBoundOn;
    else
        upBoundH = upBoundOf;
    end

    if lwBoundOn < lwBoundOf
        lwBoundH = lwBoundOn;
    else
        lwBoundH = lwBoundOf;
    end
    
    if upBoundH > upBoundM
        upBoundMH = upBoundH;
    else
        upBoundMH = upBoundM;
    end
    
    if lwBoundH < lwBoundM
        lwBoundMH = lwBoundH;
    else
        lwBoundMH = lwBoundM;
    end

    lims.audioMH = [-0.5 1.0 lwBoundMH upBoundMH];
else
    lims.audioMH = [-0.5 1.0 -100 50];
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
res.etMH   = auAn.etMH;               % 

% Raw Recorded data
res.sRate         = auAn.sRate;
res.numTrial      = auAn.numTrial;    % Total trials recorded
res.audioM        = auAn.audioM;      % Raw microphone recording (no lag correction)
res.audioH        = auAn.audioH;      % Raw headphone recording  (no lag correction)
res.allAuMHDelays = auAn.allAuMHDelays; % Vector of the delays between the NIDAQ and Audapter microphone recordings
res.allAuNiDelays = auAn.allAuNiDelays; % Vector of the delays between the NIDAQ and Audapter microphone recordings

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
res.limitsA           = lims.audioM;
res.limitsAudRes      = lims.audioAudRespMH;

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
res.limitsAMH        = lims.audioMH;

%Inflation Response
res.respVar      = auAn.respVar;
res.respVarM     = auAn.respVarM;
res.respVarSD    = auAn.respVarSD;
res.InflaStimVar = auAn.InflaStimVar;

%Some Final Output time series data
res.audioMinc = auAn.audioMinc;
res.audioHinc = auAn.audioHinc;
res.AuNiDelaysinc = auAn.AuNiDelaysinc;
end