function [auAn, auRes] = dfAnalysisAudapter(dirs, expParam, rawData, f0b, AudFlag, iRF, niAn)
% [auAn, auRes] = dfAnalysisAudapter(dirs, expParam, rawData, f0b, AudFlag, iRF, niAn)
% This function analyzes the raw audio data that was recorded by Audapter 
% in the experiments measuring changes in f0. It first does a
% pre-processing step where it identifies any experimental errors in
% production, and also identifies and corrects for any lags in recording.
% Once all the data are set up correctly and processed, they are handed off
% to a function to do the actual analysis of the audio signals. 
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
auAn.AnaType   = 'Audapter';
auAn.expType   = expParam.expType;
auAn.subject   = expParam.subject;
auAn.run       = expParam.run;
auAn.curSess   = expParam.curSess;
auAn.f0Type    = 'Praat';
auAn.f0AnaFile = [auAn.subject auAn.run 'f0Analysis' auAn.f0Type '.mat'];
auAn.gender    = expParam.gender;
auAn.AudFB     = expParam.AudFB;
auAn.AudFBSw   = expParam.AudFBSw;
auAn.bTf0b     = f0b;

fprintf('\nStarting Audapter Analysis for %s, %s with f0 of %0.2f Hz\n', auAn.subject, auAn.run, niAn.bTf0b)

% Idenitfy some Recording Variables
auAn.sRate        = expParam.sRateAnal;        % Sampling Rate of Audapter (down-sampled)
auAn.sRateNi      = niAn.sRate;                % Sampling Rate of the NIDAQ
auAn.frameLenDown = expParam.frameLen/expParam.downFact;
auAn.trialLen     = expParam.trialLen;
auAn.numSamp      = auAn.sRate*auAn.trialLen;  % Number of Samples recorded
auAn.numTrial     = expParam.numTrial;         % Number of Trials recorded
auAn.trialType    = expParam.trialType;        % Control (0), Perturbed (1)
auAn.expTrigs     = expParam.trigs(:,:,1); % Time
auAn.anaTrigs     = expParam.trigs(:,:,3); % Points (Audadapter)

auAn.time       = (0:1/auAn.sRate:(auAn.numSamp-1)/auAn.sRate)';
auAn.audioM     = [];
auAn.audioH     = [];
auAn.svIdx      = []; % Full index of trials saved against the actual trials recorded
auAn.pertIdx    = []; % The trial index pertaining to all the list of all trials saved, post-processing
auAn.contIdx    = []; % The trial index pertaining to all the list of all trials saved, post-processing
auAn.audioMSv   = [];
auAn.audioHSv   = [];
auAn.trialTypeSv = [];
auAn.expTrigsSv = [];
auAn.contTrig   = [];
auAn.pertTrig   = [];
auAn.allAuNiDelays = [];

svC = 0;
for ii = 1:auAn.numTrial
    data = rawData(ii);
    
    Mraw = data.signalIn;     % Microphone
    Hraw = data.signalOut;    % Headphones
    expTrigs = auAn.expTrigs(ii,:);
    anaTrigs = auAn.anaTrigs(ii,:);
    
    MrawNi = niAn.audioM(:,ii);   
   
    % Preprocessing step identifies time-series errors in production/recording
    [mic, head, sigDelay, preProSt] = preProcAudio(auAn, Mraw, Hraw, MrawNi, anaTrigs);
    
    auAn.audioM = cat(2, auAn.audioM, mic);
    auAn.audioH = cat(2, auAn.audioH, head);

    if preProSt.saveT == 0 %Don't save the trial :(
        fprintf('%s Trial %d not saved. %s\n', auAn.curSess, ii, preProSt.saveTmsg)
    elseif preProSt.saveT == 1 %Save the Trial
        svC = svC + 1;
        
        auAn.svIdx      = cat(1, auAn.svIdx, ii);
        auAn.expTrigsSv = cat(1, auAn.expTrigsSv, expTrigs);
        if auAn.trialType(ii) == 0
            auAn.contIdx  = cat(1, auAn.contIdx, svC);
            auAn.contTrig = cat(1, auAn.contTrig, expTrigs);
        else
            auAn.pertIdx  = cat(1, auAn.pertIdx, svC);
            auAn.pertTrig = cat(1, auAn.pertTrig, expTrigs);
        end   
        auAn.allAuNiDelays = cat(1, auAn.allAuNiDelays, sigDelay);
    end
end

%Find only the trials we care about
auAn.audioMSv      = auAn.audioM(:, auAn.svIdx);
auAn.audioHSv      = auAn.audioH(:, auAn.svIdx);
auAn.trialTypeSv   = auAn.trialType(auAn.svIdx);
auAn.numSaveTrials = length(auAn.svIdx);
auAn.numContTrials = length(auAn.contIdx);
auAn.numPertTrials = length(auAn.pertIdx);

%The Audio Analysis
f0Flag = 1;
auAn = dfAnalysisAudio(dirs, auAn, AudFlag, iRF, f0Flag);

lims  = identifyLimits(auAn);
auRes = packResults(auAn, lims);
end

function [micP, headP, AuNidelay, pp] = preProcAudio(An, micR, headR, micRNi, auTrigs)
% [micP, headP, AuNidelay, pp] = preProcAudio(An, micR, headR, micRNi, auTrigs)
% This function performs pre-processing on the time-series recorded audio 
% data before frequency analysis methods are applied. 
% 
% Inputs:
% An:      Analysis variables structure
% micR:    Raw Microphone signal (Audapter)
% headR:   Raw Headphone signal  (Audapter)
% micRNI:  Raw Microphone signal (NIDAQ)
% auTrigs: Trigger points of perturbation onset and offset (per trial; 
%
% Outputs:
% micP:      Processed Microphone signal
% headP:     Processed Headphone signal
% AuNidelay: Calculated delay between NIDAQ and Audapter
% pp:        Preprocessing result structure Reason, if any, that the trial was thrown out

micR      = double(micR);    % Convert to data type double
headR     = double(headR);   % Convert to data type double
AudFB     = An.AudFB;        % Auditory feedback type used
fs        = An.sRate;        % Sampling rate (Audapter)
fsNI      = An.sRateNi;      % Sampling rate (NIDAQ)
frameLen  = An.frameLenDown; % Frame rate of recording (After downsampling)
numSamp   = An.numSamp;      % Number of samples for length of recording

% We are going to section the audio recording from 0.5s ahead of 
% perturbation onset to 1.0s after perturbation offset.
preOn   = 0.5*fs;
postOff = 1.0*fs;

micRds     = resample(micR, fsNI, fs);
AuNidelay  = xCorrTimeLag(micRNi, micRds, fsNI); % Expected that NIDAQ will lead Audapter
AuNidelayP = AuNidelay*fs;

if strcmp(AudFB, 'Masking Noise')
    AuMHdelay = (frameLen*12)/fs;
else
    AuMHdelay = xCorrTimeLag(micR, headR, fs);   % Expected that Mic will lead Head
end
AuMHdelayP = AuMHdelay*fs;

%Adjust for delay between raw Audapter Mic and Audapter Headphones
micAuAl  = micR(1:(end-AuMHdelayP));
headAuAl = headR((AuMHdelayP+1):end); 

%Adjust for delay between Audapter and NIDAQ
micAuNi    = micAuAl(AuNidelayP:end);
headAuNi   = headAuAl(AuNidelayP:end);

%The period on either side of the pertrubation period.
audioSecSt = auTrigs(1) - preOn;
audioSecSp = auTrigs(2) + postOff;

%Section the both audio samples around section period. 
sectionInd = audioSecSt:audioSecSp;

time    = sectionInd/fs;
micSec  = micAuNi(sectionInd);
headSec = headAuNi(sectionInd);

%Find the onset of Voicing
pp = findVoiceOnsetThresh(micAuNi, fs, audioSecSt, audioSecSp);

if pp.voiceOnsetLate
    saveT    = 0;  
    saveTmsg = 'Participant started too late!!';
elseif pp.chk4Break
    saveT    = 0;
    saveTmsg = 'Participant had a voice break!!';
elseif length(micAuNi) < numSamp
    saveT    = 0;
    saveTmsg = 'Recording too short';
else
    saveT    = 1;
    saveTmsg = 'Everything is good'; 
end

timeSec = time; 
micP    = micAuNi(1:numSamp);
headP   = headAuNi(1:numSamp);

pp.saveT    = saveT;
pp.saveTmsg = saveTmsg;
end

function timeLag = xCorrTimeLag(sig1, sig2, fs)
% xCorrTimeLag(sig1, sig2, fs) is a simple function I am using to find the
% lag between two (seemingly) identical time based signals. 
% if timeLag is positive, then sig1 leads sig2. 
% if timeLag is negative, then sig1 lags sig2.

% Uses a crosscorrelation between the two signals, then finds the largest
% peak
[r, lags]    = xcorr(sig1, sig2);
[~, peakInd] = max(r);
maxLag       = lags(peakInd);
timeLag      = maxLag/fs;
timeLag      = -timeLag;
end

function pp = findVoiceOnsetThresh(audio, fs, audioSt, audioSecSp)
% findVoiceOnsetThresh is a simple script for finding the onset of voicing
% for a microphone recording. This script also makes a determination about
% whether the speaker starting speaking after a deterministic trigger
% point, or if the participant had a voice break during the production. 
%
% This returns a structure with interesting values pertaining to the
% threshold analysis, especially some notes about if the participant
% started late, or if they had a voice break.

pp.thresh   = 0.3;
pp.breakTol = 0.1;
pp.fs       = fs;
pp.audio    = audio;
pp.audioSt  = audioSt;
pp.lenSig   = length(audio);
pp.t        = 0:1/fs:(pp.lenSig-1)/fs;

[B,A] = butter(4,40/(fs/2)); %Low-pass filter under 40

%Envelope the signal removing all high frequncies. 
%This shows the general change in amplitude over time (~RMS). 
pp.env  = filter(B,A,abs(audio));  

%The largest peak in the envelope theoretically occurs during voicing
pp.maxPeak = max(pp.env); 

%I don't want to start my signal at the max value, so start lower down on
%the envelope as a threshold
pp.threshIdx     = find(pp.env > pp.thresh*pp.maxPeak); 

%The first index of the theoretical useable signal (Voice onset)
pp.voiceOnsetInd = pp.threshIdx(1);

%Check the voice onset time against when we want to start analyzing data
pp.voiceOnsetLate = audioSt < pp.voiceOnsetInd;

%The rest of the signal base the first index...are there any dead zones??
pp.fallOffLog = pp.env(pp.voiceOnsetInd:audioSecSp) < pp.thresh*pp.maxPeak;
pp.chk4Break  = sum(pp.fallOffLog) > pp.breakTol*fs; % Last longer than 300ms
end

function lims = identifyLimits(An)
% identifyLimits(An) calculates limits of analyzed data so that the limits
% are dynamic and fit the data being shown. 
%
% lims is a structure of the resultant limits to be used in plotting. 
% This function is redundant between a few different functions, and might
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
% packResults(auAn, lims) takes the results of the analysis and packages
% the important variables into a new structure that have common names
% between other analysis methods. This makes it easy to switch back and
% forth between different result structures for plotting. 

% Information about the experiment/subject
res.expType = auAn.expType;
res.subject = auAn.subject;
res.run     = auAn.run;
res.curSess = auAn.curSess;
res.AudFB   = auAn.AudFB;

res.f0Type = auAn.f0Type;
res.etMH   = auAn.etMH;

res.numTrials     = auAn.numTrial;
res.audioM        = auAn.audioM;
res.audioH        = auAn.audioH;
res.svIdx         = auAn.svIdx;
res.expTrigsSv    = auAn.expTrigsSv;
res.pertIdx       = auAn.pertIdx;     % The indices of the svIdx;
res.pertTrig      = auAn.pertTrig;
res.contIdx       = auAn.contIdx;     % The indices of the svIdx;
res.contTrig      = auAn.contTrig;
res.numSaveTrials = auAn.numSaveTrials;
res.numContTrials = auAn.numContTrials;
res.numPertTrials = auAn.numPertTrials;
res.allAuNiDelays = auAn.allAuNiDelays;

res.timef0        = auAn.timef0;
res.f0b           = auAn.f0b;

res.numContTrialsPP = auAn.numContTrialsPP;
res.numPertTrialsPP = auAn.numPertTrialsPP;
res.pertTrigPP      = auAn.pertTrigsR;

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

%NIAu Delay
res.allAuNiDelays = auAn.allAuNiDelays;
end