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
auAn.AnaType   = 'Audapter';
auAn.expType   = expParam.expType;
auAn.subject   = expParam.subject;
auAn.run       = expParam.run;
auAn.curSess   = expParam.curSess;
auAn.f0Type    = 'Praat';
auAn.f0AnaFile = [auAn.subject auAn.run 'f0Analysis.mat'];
auAn.gender    = expParam.gender;
auAn.AudFB     = expParam.AudFB;
auAn.AudFBSw   = expParam.AudFBSw;
auAn.bTf0b     = f0b;

fprintf('\nStarting Audapter Analysis for %s, %s with f0 of %0.2f Hz\n', auAn.subject, auAn.run, auAn.bTf0b)

% Idenitfy some Recording Variables
auAn.sRate        = expParam.sRateAnal;        % Sampling Rate of Audapter (down-sampled)
auAn.sRateNi      = niAn.sRate;                % Sampling Rate of the NIDAQ
auAn.frameLenDown = expParam.frameLen/expParam.downFact;
auAn.trialLen     = expParam.trialLen;         % Length of recording (s)
auAn.numSamp      = auAn.sRate*auAn.trialLen;  % Length of recording (points)

auAn.time       = (0:1/auAn.sRate:(auAn.numSamp-1)/auAn.sRate)'; % Time Vector based on numSamp
auAn.audioM     = []; % Mic data that has received temporal preprocessing (all recorded trials)
auAn.audioH     = []; % Head data that has received temporal preprocessing (all recorded trials)
auAn.numTrial   = expParam.numTrial;         % Number of trials recorded
auAn.trialType  = expParam.trialType;        % Key for identifying Control (0) & Perturbed (1) trials
auAn.expTrigs   = expParam.trigs(:,:,1); % Trigger Onset and Offset (Time) (all recorded trials)
auAn.anaTrigs   = expParam.trigs(:,:,3); % Trigger Onset and Offset (Points; Audadapter) 
auAn.removedTrialTracker = {}; % List of Trials that were thrown out during Analysis

auAn.audioMSvt     = []; % Microphone recordings for the trials saved for further analyses
auAn.audioHSvt     = []; % Headphone recordings for the trials saved for further analyses
auAn.numTrialSvt   = []; % Number of trials saved for further analyses
auAn.allIdxSvt     = []; % Vector of indicies of all recorded trials saved for further analyses.
auAn.trialTypeSvt  = []; % Key for identifying Control (0) & Perturbed (1) trials
auAn.expTrigsSvt   = []; % Trigger Onset and Offset (Time) for trials saved for further analyses
auAn.allAuNiDelays = []; % Vector of the delays between the NIDAQ and Audapter microphone recordings

auAn.numPertTrialSvt = []; % Number of perturbed trials saved for further analyses
auAn.pertIdxSvt      = []; % Vector of indicies of perturbed SAVED trials (Referencing allIdxSvt)
auAn.pertTrigSvt     = []; % Trigger Onset and Offset (Time) for perturbed SAVED trials
auAn.numContTrialSvt = []; % Number of control trials saved for further analyses
auAn.contIdxSvt      = []; % Vector of indicies of control SAVED trials (Referencing allIdxSvt)
auAn.contTrigSvt     = []; % Trigger Onset and Offset (Time) for control SAVED trials

svC = 0; % Saved Trial Count
for ii = 1:auAn.numTrial
    data = rawData(ii);       % Get the data from this trial
    
    Mraw = data.signalIn;     % Microphone
    Hraw = data.signalOut;    % Headphones
    expTrigs = auAn.expTrigs(ii,:);
    anaTrigs = auAn.anaTrigs(ii,:);
    
    if isfield(niAn, 'audioM')
        MrawNi = niAn.audioM(:,ii);          % Microphone (NIDAQ) 
    else
        MrawNi = resample(Mraw, 8000, 16000);
    end
    
    % Preprocessing step identifies time-series errors in production/recording
    [mic, head, sigDelay, preProSt] = preProcAudio(auAn, Mraw, Hraw, MrawNi, anaTrigs);
    
    auAn.audioM = cat(2, auAn.audioM, mic);  % Save all trials, regardless of eventual analysis
    auAn.audioH = cat(2, auAn.audioH, head); % Save all trials, regardless of eventual analysis

    if preProSt.saveT == 0     % Don't save the trial :(
        fprintf('%s Trial %d not saved. %s\n', auAn.curSess, ii, preProSt.saveTmsg)
        removedTrial = {['Trial ' num2str(ii)], preProSt.saveTmsg};
        auAn.removedTrialTracker = cat(1, auAn.removedTrialTracker, removedTrial);
    elseif preProSt.saveT == 1 % Save the Trial
        svC = svC + 1; % Iterate Saved Trial Count
        
        auAn.allIdxSvt   = cat(1, auAn.allIdxSvt, ii); % Save the experimental index
        auAn.expTrigsSvt = cat(1, auAn.expTrigsSvt, expTrigs); % Save the triggers from this index
        if auAn.trialType(ii) == 0
            auAn.contIdxSvt  = cat(1, auAn.contIdxSvt, svC); %
            auAn.contTrigSvt = cat(1, auAn.contTrigSvt, expTrigs); %
        else
            auAn.pertIdxSvt  = cat(1, auAn.pertIdxSvt, svC); %
            auAn.pertTrigSvt = cat(1, auAn.pertTrigSvt, expTrigs); %
        end   
        auAn.allAuNiDelays = cat(1, auAn.allAuNiDelays, sigDelay);
    end
end

% Find only the trials we care about
auAn.audioMSvt     = auAn.audioM(:, auAn.allIdxSvt); % Grabbing the recorded audio based on the saved indices
auAn.audioHSvt     = auAn.audioH(:, auAn.allIdxSvt); % Grabbing the recorded audio based on the saved indices
auAn.trialTypeSvt  = auAn.trialType(auAn.allIdxSvt); % The order of trial type based on the saved trials post-tP

auAn.numTrialSvt     = length(auAn.allIdxSvt);
auAn.numPertTrialSvt = length(auAn.pertIdxSvt);
auAn.numContTrialSvt = length(auAn.contIdxSvt);

% The Audio Analysis
f0Flag = 0;
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
if AuNidelayP > 0 % As long as the delay is non 0
    micAuNi    = micAuAl(AuNidelayP:end);
    headAuNi   = headAuAl(AuNidelayP:end);
else
    micAuNi    = micAuAl;
    headAuNi   = headAuAl;
end

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
res.expType = auAn.expType;
res.subject = auAn.subject;
res.run     = auAn.run;
res.curSess = auAn.curSess;
res.AudFB   = auAn.AudFB;

res.f0Type = auAn.f0Type;
res.etMH   = auAn.etMH;

% Raw Recorded data
res.numTrial      = auAn.numTrial;    % Total trials recorded
res.audioM        = auAn.audioM;      % Raw microphone recording (no lag correction)
res.audioH        = auAn.audioH;      % Raw headphone recording  (no lag correction)
res.removedTrialTracker = auAn.removedTrialTracker;

% Post temporal processing data
res.numTrialSvt   = auAn.numTrialSvt;   % Number of trials saved (post temporal processing)
res.allIdxSvt     = auAn.allIdxSvt;     % Vector of indicies of recorded trials saved (post temporal processing)
res.trialTypeSvt  = auAn.trialTypeSvt;  % Key for identifying Control (0) & Perturbed (1) trials (post temporal processing)
res.expTrigsSvt   = auAn.expTrigsSvt;   % Trigger Onset and Offset (Time) for trials saved (post temporal processing)
res.allAuNiDelays = auAn.allAuNiDelays; % Vector of the delays between the NIDAQ and Audapter microphone recordings

res.numPertTrialSvt = auAn.numPertTrialSvt; % Number of perturbed trials saved (post temporal processing)
res.pertIdxSvt      = auAn.pertIdxSvt;      % Vector of indicies of perturbed trials (Referencing allIdxSvt) (post temporal processing)
res.pertTrigSvt     = auAn.pertTrigSvt;     % Trigger Onset and Offset (Time) for perturbed trials (post temporal processing)
res.numContTrialSvt = auAn.numContTrialSvt; % Number of control trials saved (post temporal processing)
res.contIdxSvt      = auAn.contIdxSvt;      % Vector of indicies of control trials (Referencing allIdxSvt) (post temporal processing)
res.contTrigSvt     = auAn.contTrigSvt;     % Trigger Onset and Offset (Time) for control trials (post temporal processing)

% Audio f0 analysis
res.timef0        = auAn.timef0;
res.f0b           = auAn.f0b;

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

%NIAu Delay
res.allAuNiDelays = auAn.allAuNiDelays;
end