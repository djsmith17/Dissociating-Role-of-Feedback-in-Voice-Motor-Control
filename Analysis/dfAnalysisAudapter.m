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

svC = 0; % Saved Trial Count
for ii = 1:auAn.numTrial
    data = rawData(ii);       % Get the data from this trial
    
    Mraw = data.signalIn;     % Microphone
    Hraw = data.signalOut;    % Headphones
    expTrigs = auAn.expTrigs(ii,:);
    anaTrigs = auAn.anaTrigs(ii,:);
    
    MrawNi = niAn.audioM(:,ii); % Microphone (NIDAQ) 
   
    % Preprocessing step identifies time-series errors in production/recording
    [mic, head, sigDelay, preProSt] = preProcAudio(auAn, Mraw, Hraw, MrawNi, anaTrigs);
    
    auAn.audioM = cat(2, auAn.audioM, mic);  % Save all trials, regardless of eventual analysis
    auAn.audioH = cat(2, auAn.audioH, head); % Save all trials, regardless of eventual analysis

    if preProSt.saveT == 0     % Don't save the trial :(
        fprintf('%s Trial %d not saved. %s\n', auAn.curSess, ii, preProSt.saveTmsg)
    elseif preProSt.saveT == 1 % Save the Trial
        svC = svC + 1; % Iterate Saved Trial Count
        
        auAn.svIdx      = cat(1, auAn.svIdx, ii); % Save the experimental index (Which experimental trial?)
        auAn.expTrigsSv = cat(1, auAn.expTrigsSv, expTrigs); % Save the triggers from this index
        if auAn.trialType(ii) == 0
            auAn.contIdx  = cat(1, auAn.contIdx, svC); %
            auAn.contTrig = cat(1, auAn.contTrig, expTrigs); %
        else
            auAn.pertIdx  = cat(1, auAn.pertIdx, svC); %
            auAn.pertTrig = cat(1, auAn.pertTrig, expTrigs); %
        end   
        auAn.allAuNiDelays = cat(1, auAn.allAuNiDelays, sigDelay);
    end
end

% Find only the trials we care about
auAn.audioMSv      = auAn.audioM(:, auAn.svIdx); 
auAn.audioHSv      = auAn.audioH(:, auAn.svIdx);
auAn.trialTypeSv   = auAn.trialType(auAn.svIdx); % The order of trial type based on the saved trials post-tP
auAn.numSaveTrials = length(auAn.svIdx);       % # of total trials saved after temporal processing
auAn.numContTrials = length(auAn.contIdx);     % # of control trials saved after temporal processing
auAn.numPertTrials = length(auAn.pertIdx);     % # of perturbed trial saved after temporal processing

% The Audio Analysis
f0Flag = 1;
auAn = dfAnalysisAudio(dirs, auAn, AudFlag, iRF, f0Flag);

lims  = identifyLimits(auAn);
auRes = packResults(auAn, lims);
end

function [micP, headP, AuNidelay, pp] = preProcAudio(An, micR, headR, micRNi, auTrigs)
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
% AuNidelay: Calculated delay between NIDAQ and Audapter
% pp:        Preprocessing results structure. This has information
%            regarding the envelope of the recorded audio file, and 
%            if the participant started late, or has a voice break. 

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
AuNidelay  = xCorrTimeLag(micRNi, micRds, fsNI); % Expect NIDAQ leads Audapter
AuNidelayP = AuNidelay*fs;

if strcmp(AudFB, 'Masking Noise')
    AuMHdelay = (frameLen*12)/fs;
else
    AuMHdelay = xCorrTimeLag(micR, headR, fs);   % Expect Mic leads Head
end
AuMHdelayP = AuMHdelay*fs;

% Adjust for delay between raw Audapter Mic and Audapter Headphones
micAuAl  = micR(1:(end-AuMHdelayP));
headAuAl = headR((AuMHdelayP+1):end); 

% Adjust for delay between Audapter and NIDAQ
micAuNi    = micAuAl(AuNidelayP:end);
headAuNi   = headAuAl(AuNidelayP:end);

% Audio points on either side of the perturbation period.
audioSecSt = auTrigs(1) - preOn;
audioSecSp = auTrigs(2) + postOff;

% Find the onset of voicing
pp = findVoiceOnset(micAuNi, fs, audioSecSt, audioSecSp);

if pp.voiceOnsetLate
    saveT    = 0;  
    saveTmsg = 'Participant started too late!!';
elseif pp.chk4Break
    saveT    = 0;
    saveTmsg = 'Participant had a voice break!!';
elseif length(micAuNi) < numSamp
    saveT    = 0;
    saveTmsg = 'Recording too short!!';
else
    saveT    = 1;
    saveTmsg = 'Everything is good'; 
end

% Grab the full numSamp so they can be concatenated cleanly
micP    = micAuNi(1:numSamp);
headP   = headAuNi(1:numSamp);

pp.saveT    = saveT;    % Save trial or no?
pp.saveTmsg = saveTmsg; % Reason, if any the trial was thrown out
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

function pp = findVoiceOnset(audio, fs, audioSt, audioSp)
% pp = findVoiceOnset(audio, fs, audioSt, audioSecSp) identifies
% onset of voice by the envelope of a microphone recording (audio).
% Based on a specific point (audioSt), this script identifies if the 
% participant started production late. Then the script identifies any
% points where the participant might have had a voice break, up to the end
% of interesting data will be collected (audioSp). 
%
% This returns a structure (pp) with the results of the voice onset 
% detection methods, and boolean values for whether the participant started
% late (pp.voiceOnsetLate) or if they had a voice break (pp.chk4Break)

pp.thresh   = 0.3; % Threshold of Decimal amount of full peak height
pp.breakTol = 0.1; % Voice Break Tolerance; Time (s)
pp.envCutOf = 40;  % Cutoff frequency for enveloping the audio
pp.fs       = fs;
pp.audio    = audio;
pp.audioSt  = audioSt;
pp.audioSp  = audioSp;
pp.lenSig   = length(audio);
pp.t        = 0:1/fs:(pp.lenSig-1)/fs; % Not used, but useful if debugging

% 4th order low-pass butter filter settings
[B,A] = butter(4, pp.envCutOf/(fs/2));

% Envelope the signal by low-pass filtering (change in amplitude/time ~RMS)
pp.env  = filter(B,A,abs(audio));  

% Largest peak in the envelope theoretically occurs during voicing
pp.maxPeak = max(pp.env); 

% Find values that are within threshold of max 'voicing' value
pp.threshIdx     = find(pp.env > pp.thresh*pp.maxPeak); 

% First index of the theoretical useable signal (Voice onset)
pp.voiceOnsetInd = pp.threshIdx(1);

% Check the voice onset time against when we want to start analyzing data
pp.voiceOnsetLate = audioSt < pp.voiceOnsetInd;

% Check the whole The rest of the signal base the first index...are there any dead zones??
pp.fallOffLog = pp.env(pp.voiceOnsetInd:audioSp) < pp.thresh*pp.maxPeak;
pp.chk4Break  = sum(pp.fallOffLog) > pp.breakTol*fs; % Last longer than break tolerance
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