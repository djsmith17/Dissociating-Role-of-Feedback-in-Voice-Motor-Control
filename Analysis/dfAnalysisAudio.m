function An = dfAnalysisAudio(dirs, An, AudFlag, varargin)
% An = dfAnalysisAudio(dirs, An, AudFlag, varargin) performs frequency 
% analyses on recorded audio data to measure the pitch contour on a number
% of recorded trials. The actual f0 measurement is done in the sub-function
% signalFrequencyAnalysis. The remainder of this script processes the
% f0-trace (for both Mic/Head) with the following steps:
% 1: Identifying the baseline f0 value for each trial
% 2: Normalizing each f0-trace by its baseline f0 value
% 3: Parsing trials between perturbed and control trial types
% 4: Sectioning each trial around the trigger values
% 5: Taking the mean section of trials of a type
% 6: Performing (if applicable) audio dynamics of the mean f0-trace
%
% INPUTS
% dirs:    Structure of the file directories set for the computer in use.
% An:      Structure of analysis variables to be shared between analysis
%          steps.
% AudFlag: Boolean flag for running frequency analysis
% adFlag:  Audio Dynamics Flag. Analyze changes in f0 following triggers
% f0Flag:  Perform f0-trace calculation, or load previous version?
%
% OUTPUTS
% An:      Structure of analysis variables to be shared between analysis
%          steps.
% 
% This function has the following subfunctions
% -initAudVar
% -setFreqAnalVar
% -signalFrequencyAnalysis
% -identifyf0Bounds
% -SimpleAutoCorr
% -smoothf0
% -normf0
% -parseTrialTypes
% -sectionData
% -round2match
% -meanAudioData
% -InflationResponse
% -InflationResponseStruct
% -PitchShiftReflexResponse
%
% This function expects that the input An contains the following fields
% -An.curSess
% -An.f0Type
% -An.f0AnaFile
% -An.f0b
% -An.gender
% -An.sRate
% -An.allIdxPreProc
% -An.audioM
% -An.audioH
% -An.expTrigs
% -An.trialType
%
% Author: Dante J Smith
% Updated: 04/06/2022

if isempty(varargin)
    adFlag = 0; 
    f0Flag = 0;
elseif length(varargin) == 1
    adFlag = varargin{1};
    f0Flag = 0;
else
    adFlag = varargin{1};
    f0Flag = varargin{2};
end

% Initalize analysis variables
An = initAudVar(An);

if AudFlag == 1
    fprintf('\nStarting Pitch Analysis\n')

    An.audioMSvt    = An.audioM(:, An.allIdxPreProc);
    An.audioHSvt    = An.audioH(:, An.allIdxPreProc);
    An.expTrigsSvt  = An.expTrigs(An.allIdxPreProc, :);
    An.trialTypeSvt = An.trialType(An.allIdxPreProc);
    An.numTrialSvt  = length(An.allIdxPreProc);
 
    freqVar.numSamp = length(An.audioMSvt);
    freqVar.sRate   = An.sRate;
    freqVar.f0b     = An.f0b;
    freqVar.gender  = An.gender;
    
    % File where to save/find f0 trace analysis
    dirs.audiof0AnalysisFile = fullfile(dirs.SavResultsDir, An.f0AnaFile);
    
    % Select method of f0 trace calculation
    if strcmp(An.f0Type, 'Praat') == 1
        freqVar.f0AnalysisType = 1;
    else
        freqVar.f0AnalysisType = 2;
    end
        
    % f0-trace measurement from time-series data 
    if exist(dirs.audiof0AnalysisFile, 'file') == 0 || f0Flag == 1
        % Perform signal frequency analysis on audio signals = f0 traces
        [f0A.timef0, f0A.audioMf0, f0A.expTrigsf0, f0A.etM, f0A.fV] = signalFrequencyAnalysis(dirs, freqVar, An.audioMSvt, An.expTrigsSvt);
        [f0A.timef0, f0A.audioHf0, f0A.expTrigsf0, f0A.etH, f0A.fV] = signalFrequencyAnalysis(dirs, freqVar, An.audioHSvt, An.expTrigsSvt);        
        save(dirs.audiof0AnalysisFile, 'f0A')
    else
        % Load previously analyzed f0 traces
        load(dirs.audiof0AnalysisFile)
    end

    An.timef0     = f0A.timef0;
    An.expTrigsf0 = f0A.expTrigsf0;
    An.audioMf0   = f0A.audioMf0;
    An.audioHf0   = f0A.audioHf0;
    An.etMH       = f0A.etM + f0A.etH; % Minutes
    An.fV         = f0A.fV;
    
    % Set up time information for sectioned-data
    preEveT = -0.5;
    posEveT = 1.0;
    eveTLen = posEveT - preEveT;
    numSampSec = eveTLen/An.fV.win + 1; % Number of samples in sectioned trial
    
    % Time vector corresponding to the sectioned signals
    An.secTime = linspace(preEveT, posEveT, numSampSec);
    
    % Smooth the f0 data
%     An.audioMf0S   = smoothf0(An.audioMf0);
%     An.audioHf0S   = smoothf0(An.audioHf0);

    % Section Audio with all trials...before parsing, and post-processing
    An.audioMf0SecAll = sectionData(An.timef0, An.audioMf0, An.expTrigsf0, numSampSec);
    An.audioHf0SecAll = sectionData(An.timef0, An.audioHf0, An.expTrigsf0, numSampSec);
   
    % Find the value of f0 during the perPert period for each trial
    prePert      = (An.secTime <= 0); % SecTime is aligned for SecTime = 0 to be Onset of pert
    An.trialf0   = nanmean(An.audioMf0SecAll(prePert,:,1),1); % Per-trial baseline f0   
    An.trialf0M  = nanmean(An.trialf0);                       % Mean trial baseline f0
    
    % Normalize f0 traces by trial f0b and convert to cents
    An.audioMf0_norm = normf0(An.audioMf0, An.trialf0);
    An.audioHf0_norm = normf0(An.audioHf0, An.trialf0);
    
    % Identify trials (mic) where participant had vocal fry
    % aka: f0 pitch miscalculation
    for ii = 1:An.numTrialSvt
        
        expTrig = An.expTrigsf0(ii, :);
        svIdc   = An.allIdxPreProc(ii);
        
        preSt = expTrig(1) - 0.5; % Pre-Start
        posSp = expTrig(2) + 1.0; % Post-Stop
        
        % Only consider points around where the action happens
        % Ignore the fringes of the recording
        timeInd = (An.timef0 > preSt & An.timef0 < posSp);
        mic     = An.audioMf0_norm(timeInd, ii);

        % Any points where pitch is > 500 // < -500 cents?
        MisCalcf0 = find(mic >= 500 | mic <=  -500);

        if ~isempty(MisCalcf0)
            if An.trialTypeSvt(ii) == 0
                type = 'Control';
            else
                type = 'Perturbed';
            end       
            fprintf('%s Trial %d (%s) excluded due to Miscalculated Pitch Trace\n', An.curSess, svIdc, type)

            removedTrial = {['Trial ' num2str(svIdc)], 'Miscalculated pitch Trace'};
            An.removedTrialTracker = cat(1, An.removedTrialTracker, removedTrial);
        else
            An.subSvIdx = cat(1, An.subSvIdx, ii); % Keep the index of the trials that are not removed
        end
    end
%%%%%%%%    
%     fs = 1/An.fV.win;
%     trialVar.coder = 'Pipe';
%     trialVar.curSess = 'Pipe';
%     trialVar.sigType = 'f0';
%     trialVar.units   = 'cents';
%     trialVar.itrType = '';
%     secf0 = dfSectionDataOrg(An.timef0, An.audioMf0, An.expTrigsf0, fs, trialVar);
%     
%     % Section raw f0 around onset and offset
%     secf0.sigsSec  = secf0.sectionData(secf0.sigs);
%     
%     % Identify baseline values
%     secf0 = secf0.identifyBaselineValues(secf0.sigsSec);
%     
%     % Convert to cents
%     secf0 = secf0.convertCentsData();
%     
%     % Section converted f0 around onset and offset
%     secf0.sigsNormSec = secf0.sectionData(secf0.sigsNorm);
%     
%     % Quality check trial
%     secf0 = secf0.qualityCheckData(secf0.sigsNormSec, An.allIdxPreProc);
%     
%     % Separate saved trials
%     secf0.sigsNormSecSv = secf0.sigsNormSec(:,secf0.svIdx,:);
%     pertTrialsIdx = An.trialTypeSvt(secf0.svIdx) == 1;
%     secf0P = secf0;
%     secf0P.sigsNormSecSv = secf0.sigsNormSecSv(:, pertTrialsIdx, :);
%     
%     % Mean the trials
%     secf0P.sigsSecM = secf0P.meanData(secf0P.sigsNormSecSv);
%     
%     % Identify Inflation Response
%     secf0P.sigsDynamics = secf0P.InflationResponse;
%     
%     % Identify limits of the mean trials
%     secf0P = secf0P.identifyBounds;
%     
%     % Draw the Onset Figure
%     secf0P.drawSigsSecM_Onset(2);
    
%%%%%%%%    
    % Parse out the trials we are saving (not removed)
    An.svf0Idx       = An.allIdxPreProc(An.subSvIdx);
    An.expTrigsf0Sv  = An.expTrigsf0(An.subSvIdx, :);
    An.trialTypef0Sv = An.trialTypeSvt(An.subSvIdx);
    An.contf0Idx     = An.trialTypef0Sv == 0;
    An.pertf0Idx     = An.trialTypef0Sv == 1;
    
    An.audioMf0sv      = An.audioMf0_norm(:, An.subSvIdx);
    An.audioHf0sv      = An.audioHf0_norm(:, An.subSvIdx); 
    An.numTrialsPP     = length(An.subSvIdx);
    An.numPertTrialsPP = sum(An.pertf0Idx);
    An.numContTrialsPP = sum(An.contf0Idx);

    %Parse between Perturbed and Control trials
    An.pertTrigsR  = An.expTrigsf0Sv(An.pertf0Idx,:);
    An.audioMf0p   = parseTrialTypes(An.audioMf0sv, An.pertf0Idx);
    An.audioHf0p   = parseTrialTypes(An.audioHf0sv, An.pertf0Idx);
    An.contTrigsR  = An.expTrigsf0Sv(An.contf0Idx,:);
    An.audioMf0c   = parseTrialTypes(An.audioMf0sv, An.contf0Idx);
    An.audioHf0c   = parseTrialTypes(An.audioHf0sv, An.contf0Idx);
    
    %Section the data around onset and offset
    An.audioMf0_Secp = sectionData(An.timef0, An.audioMf0p, An.pertTrigsR, numSampSec);
    An.audioHf0_Secp = sectionData(An.timef0, An.audioHf0p, An.pertTrigsR, numSampSec);
    An.audioMf0_Secc = sectionData(An.timef0, An.audioMf0c, An.contTrigsR, numSampSec);
    An.audioHf0_Secc = sectionData(An.timef0, An.audioHf0c, An.contTrigsR, numSampSec);

    %Mean around the onset and offset
    An.audioMf0_meanp = meanAudioData(An.audioMf0_Secp);
    An.audioHf0_meanp = meanAudioData(An.audioHf0_Secp);
    An.audioMf0_meanc = meanAudioData(An.audioMf0_Secc);
    An.audioHf0_meanc = meanAudioData(An.audioHf0_Secc); 
    
    if adFlag == 1     % Set in RunSubjAnalysis
        An.audioDynamics = InflationResponse(An.secTime, An.audioMf0_meanp); % Audio Response to Somatosensory Pert
    elseif adFlag == 2 % Set in RunSubjAnalysis
        An.audioDynamics = PitchShiftReflexResponse(An.secTime, An.audioMf0_meanp); % Audio Response to Auditory Pert
    end
end
end

function An = initAudVar(An)
% An = initAudVar(An) initialize variables used in the analysis of 
% recorded audio (Mic and Head).
%
% INPUTS
% An: Structure of analysis variables to be shared between analysis
%     steps.
%
% OUTPUTS
% An: Structure of analysis variables to be shared between analysis
%     steps. Updated with initalized values.

An.timef0         = []; % Time vector of audio samples recorded
An.audioMf0       = []; % Vector of raw Microphone audio
An.audioHf0       = []; % Vector of raw Headphone audio
An.expTrigsf0     = []; % ExpTrigs that have been slightly time shifted to match sampling rate of f0 analysis
An.etMH           = []; % Float value of the amount of elapsed time to perform analysis of Microphone and headphone data
An.fV             = []; % Structure of frequency analysis variables

An.audioMf0S      = [];
An.audioHf0S      = [];

An.secTime        = [];
An.audioMf0SecAll = [];
An.audioHf0SecAll = [];

An.trialf0        = []; % Per-Trial f0 from the period pre-perturbation
An.trialf0M       = []; % Mean-Trial f0 from the period per-perturbation
An.audioMf0_norm  = []; % Normalized mic data (cents)
An.audioHf0_norm  = []; % Normalized head data (cents)

An.subSvIdx       = [];
An.svf0Idx        = [];
An.expTrigsf0Sv   = [];
An.trialTypef0Sv  = [];
An.pertf0Idx      = [];
An.contf0Idx      = [];

An.audioMf0sv      = [];
An.audioHf0sv      = []; 
An.numTrialsPP     = [];
An.numPertTrialsPP = [];
An.numContTrialsPP = [];

An.pertTrigsR    = [];
An.contTrigsR    = [];
An.audioMf0p     = []; % Perturbed mic data (norm)
An.audioHf0p     = []; % Pertrubed head data (norm)
An.audioMf0c     = []; % Control mic data (norm)
An.audioHf0c     = []; % Control head data (norm)

An.audioMf0_Secp  = []; 
An.audioHf0_Secp  = [];
An.audioMf0_Secc  = []; 
An.audioHf0_Secc  = [];

An.audioMf0_meanp = []; 
An.audioHf0_meanp = [];
An.audioMf0_meanc = []; 
An.audioHf0_meanc = [];

An.audioDynamics  = []; % Structure of Audio Dynamics results
end

function fV = setFreqAnalVar(freqVar)
% fV = setFreqAnalVar(freqVar) sets frequency analysis variables for
% analysis of audio signals.
%
% INPUTS
% freqVar: Structure of analysis variables to be used when calculating fV
%
% OUTPUTS
% fV: Structure of frequency analysis variables to be used in frequency
%     analysis of audio signals

% Pull out basic recording 
fV.anaFlag    = freqVar.f0AnalysisType;
fV.sRate      = freqVar.sRate;
fV.numSamp    = freqVar.numSamp;

% Set analysis variables
fV.f0b        = freqVar.f0b;
fV.gender     = freqVar.gender;
fV.f0Bounds   = identifyf0Bounds(fV.gender, fV.f0b);

fV.time       = (0:1/fV.sRate:(fV.numSamp-1)/fV.sRate)'; % Time vector for full mic

fV.freqCutOff = 400;        % Cut Off Frequency.
fV.win        = 0.001;      % Sampling window length (time; s)
fV.fsA        = 1/fV.win;   % 
fV.winP       = fV.win*fV.sRate; % Sampling window length (points; n)
fV.pOV        = 0.00;       % 0% overlap
fV.tStepP     = round(fV.winP*(1-fV.pOV)); %Number of windows steps over the sampling period

fV.trialWin = round(1:fV.tStepP:(fV.numSamp-fV.winP)); % Window start frames based on length of voice onset mic
fV.numWin   = length(fV.trialWin);                     % Number of windows based on WinSt

fV.roundFact = fV.sRate/fV.tStepP; % Rounding factor for frequency rounding
fV.winHalf   = fV.win/2;           % Half a sampling window length (s)
end

function [timef0, audiof0, expTrigsR, elapsed_Time, fV] = signalFrequencyAnalysis(dirs, freqVar, audio, expTrig)
%[timef0, audiof0, expTrigsR, elapsed_Time, fV] = signalFrequencyAnalysis
% (dirs, freqVar, audio, expTrig) performs the frequency analysis required
% for these expertiments. Specifically this function offers two methods two
% pull out the fundamental frequency (f0) from the recorded voice
% recordings. The first method is to use Praat methods to pull out the f0
% from the voice recordings and saving them in separate file. The second
% method is to perform an autocorrelation on the audio signal to pull out
% the f0 trace. 
% 
% INPUTS
% dirs:         Structure of the file directories set for the computer in use
% freqVar:      Structure of variables used to set the frequency variables
% audio:        Matrix of audio signals to perform analysis on
% expTrig:      Matrix of experimental triggers
%
% OUTPUTS
% timef0:       Vector of time points corresponding to the f0-traces
% audiof0:      Matrix of fundamental frequency measures taken from the
%               time-series audio data.
% expTrigsR:    Vector of experimental triggers rounded
% elapsed_Time: Amount of time that has elapsed during this analysis step
% fV:           Structure of frequency analysis variables to be used in 
%               frequency analysis of audio signals


elapsed_Time_Start = tic;
[~, numTrial] = size(audio);

fV = setFreqAnalVar(freqVar);

if fV.anaFlag == 1    
    [timef0, audiof0] = dfCalcf0Praat(dirs, audio, fV.sRate, fV.win, fV.f0Bounds);
else
    audiof0 = [];
    for j = 1:numTrial % Trial by Trial             
        time  = fV.time;
        voice = audio(:,j);
         
        timef0  = []; 
        voicef0 = [];
        for i = 1:fV.numWin
            winIdx  = fV.trialWin(i):fV.trialWin(i)+ fV.winP-1;
            
            timeWin  = round(mean(time(winIdx)),4);
            voiceWin = voice(winIdx);
            
            % What is the f0??
%             f0Win = dfCalcf0Chile(voiceWin, fs);
%             if isnan(f0Win); f0Win = bTf0b; end
            
            f0Win = simpleAutoCorr(voiceWin, fV);
            
            timef0  = cat(1, timef0, timeWin);
            voicef0 = cat(1, voicef0, f0Win);
        end    
        audiof0 = cat(2, audiof0, voicef0);
    end
end

expTrigsR = round2matchfs(expTrig, fV.roundFact, fV.winHalf);

elapsed_Time = toc(elapsed_Time_Start)/60;
% fprintf('f0 Analysis Elapsed Time: %f (min)\n', elapsed_Time)
end

function bounds = identifyf0Bounds(gender, f0b)
% bounds = identifyf0Bounds(gender, f0b) is a simple script to set the
% upper and lower bounds of fundamental frequency (f0) for a given
% pariticipant based on their birth sex (gender in short form). This is the
% possible range that a given participant's f0 may fluctuate during a given
% experimental task. These bounds will be used to selectively tune the
% Audapter algorithms which perturb auditory feedback. The default bounds
% are set based on literature review and adjusted in cases of extreme
% variation in the recorded baseline fundamental frequency of the
% participant,
%
% INPUTS
% gender: String value of birth sex ('Male' or 'Female') recorded from 
% baseline assessments
% f0b: Integer value of fundamental frequency of the participant recorded 
% during baseline assessments
%
% OUTPUTS
% bounds: Vector of lower(1) and upper(2) limits of fundamental frequency
% (Hz) to be used for auditory feedback perturbations for the selected
% participant.

% Default fundamental frequency bounds identified from literature review
defaultMale   = [75 300];  % (Hz)
defaultFemale = [100 500]; % (Hz)

switch gender
    case 'male'
        if (f0b/2) < defaultMale(1) % Especially low-pitch Male
            bounds = [25 250]; % (Hz)
        else
            bounds = defaultMale;
        end
        
    case 'female'
        if (f0b*2) > defaultFemale(2) % Especially high-pitch Female
            bounds = [200 600]; % (Hz)
        else
            bounds = defaultFemale;
        end
end
end

function f0Win = simpleAutoCorr(voice, fV)
% f0Win = simpleAutoCorr(voice, fV) is a very simple autocorrelation method
% for finding the fundamental frequency during the frequency analysis.
%
% INPUTS
% voice: vector of audio/voice signal
% fV: Structure of frequency variables
%
% OUTPUTS
% f0Win: Frequency of the highest peak in the recording

fs            = fV.sRate;               % Sampling Rate
win           = fV.winP;                % Sampling Window
[autoC, lags] = autocorr(voice, win-1); % Perform autocorrelation
[pks, pkInd]  = findpeaks(autoC);       % findoeaks from autocorrelation

if isempty(pks)
    FLag      = lags(end);
else
    [~, mxInd]= max(pks);
    FLag      = lags(pkInd(mxInd));
end
per           = FLag/fs;
f0Win         = 1/per;
end

function audioS = smoothf0(audio)
% To be deleted
[~, numTrial] = size(audio);

audioS = [];
for ii = 1:numTrial
    audioSmooth = smooth(audio(:,ii), 50);   % 10 sample length smoothing
    audioS      = cat(2, audioS, audioSmooth);
end
end

function audio_norm = normf0(audio, trialf0)
% audio_norm = normf0(audio, trialf0) takes a matrix of audio signals 
% (audio) of size numSamp x numTrial and normalizes each trial by the 
% baseline f0 caluclated for that trial, stored in the vector trialf0. 
%
% INPUTS
% audio: matrix of audio signals (numSamp x numTrial)
% trialf0: vector of baseline f0 values (numTrial x 1)
%
% OUTPUTS
% audio_norm: matrix of normalized audio signals (numSamp x numTrial)

[~, numTrial] = size(audio);

audio_norm = [];
for ii = 1:numTrial
    audio_trial = 1200*log2(audio(:,ii)./trialf0(ii));
    audio_norm  = cat(2, audio_norm, audio_trial);
end
end

function signalParse = parseTrialTypes(signal, idx)
% signalParse = parseTrialTypes(signal, idx) parses individual trials out 
% of a large matrix of recordings of size numSamp x numTrial. 
% (idx) is a vector of the indices to parse out.
% Why did you make this a function? Get over it.
%
% INPUTS
% signal: matrix of audio signals (numSamp x numTrial)
% idx: Index to sort out
%
% OUTPUTS
% signalParse: Matrix of parsed signals (numSamp x trials of a type)

signalParse = signal(:, idx);
end

function secSigs = sectionData(time, sigs, trigs, numSampSec)
% secSigs = sectionData(time, sigs, trigs, numSampSec) sections time-series
% data around important points in time.
% 
% INPUTS
% time:  Vector of time points (numSamp)
% sigs:  Matrix of signals to be sectioned (numSamp x numTrial)
% trigs: Onset and Offset time tiggers (numTrial x 2)
% numSampSec: Number of Samples in sectioned-data
%
% OUTPUTS
% secSigs: 3D mat of sectioned sigs (numSampSec x numTrial x event)
%          The 1st 3D layer are Onset Sections
%          The 2nd 3D later are Offset Sections

[~, numTrial] = size(sigs);
preEveT  = -0.5;

% Sectioned-data vector of points
pVec = linspace(0, numSampSec-1, numSampSec);

OnsetSecs  = [];
OffsetSecs = [];
if numTrial > 0
    for ii = 1:numTrial
        thisSig = sigs(:, ii);
        thisSig = padarray(thisSig, 16000, 0, 'post');
        
        OnsetT   = trigs(ii, 1); % Onset time
        OffsetT  = trigs(ii, 2); % Offset time

        OnsetTSt = round(OnsetT + preEveT, 3);   % PreOnset time, rounded to nearest ms
        OnsetTStLeast = find(time <= OnsetTSt);
        OnsetSpan = OnsetTStLeast(end) + pVec; % Indices corresponding to Onset period

        OffsetTSt = round(OffsetT + preEveT, 3); % PreOffset time, rounded to nearest ms
        OffsetTStLeast = find(time <= OffsetTSt);
        OffsetSpan = OffsetTStLeast(end) + pVec; % Indices corresponding to Offset period

        OnsetSec  = thisSig(OnsetSpan);  % Data sectioned around Onset
        OffsetSec = thisSig(OffsetSpan); % Data sectioned around Offset

        OnsetSecs  = cat(2, OnsetSecs, OnsetSec);   % Sectioned signal onsets concatenated
        OffsetSecs = cat(2, OffsetSecs, OffsetSec); % Sectioned signal offsets concatenated
    end
end

secSigs(:,:,1) = OnsetSecs;  % 1st 3D layer
secSigs(:,:,2) = OffsetSecs; % 2nd 3D layer
end

function y = round2matchfs(x, rFact, winHalf)
% y = round2matchfs(x, rFact, winHalf) rounds the value of x to a value y,
% which is the closet number to x that matches the time step for a given
% analysis windowing methods.
%
% INPUTS
% x: Decimal value of frequency. X can be a single number, a vector of 
%          numbers or a matrix of numbers.
% rfact: sampling rate divided by the sample size
% winHalf: Half a sampling window (s)
%
% OUTPUTS
% y: Integer value of frequency rounded from input x

y = round((x-winHalf).*rFact)./rFact + winHalf;
end

function meanAudio = meanAudioData(secAudio)
% meanAudio = meanAudioData(secAudio) performs the averaging and 
% Standard Error measurements for the sectioned audio trials.
%
% INPUTS
% secAudio: 3D matrix of sectioned sigs (numSampSec x numTrial x event)
%           The 1st 3D layer are Onset Sections
%           The 2nd 3D later are Offset Sections
% OUTPUTS
% meanAudio: 2D Matrix of mean section sigs (numSampSec x numColumns)
%            The columns included in this output matrix are the following:
%            meanAudio(:,1) = mean Onset pitch contour
%            meanAudio(:,2) = Standard Error of the mean Onset Pitch Contour
%            meanAudio(:,3) = mean Offset pitch contour
%            meanAudio(:,4) = Standard Error of the mean Offset Pitch Contour

OnsetSecs  = secAudio(:,:,1);
OffsetSecs = secAudio(:,:,2);
[~, numTrial] = size(OnsetSecs);

meanOnset  = nanmean(OnsetSecs, 2);  % across columns
meanOffset = nanmean(OffsetSecs, 2); % across columns

stdOnset   = nanstd(OnsetSecs, 0, 2);  % across columns
stdOffset  = nanstd(OffsetSecs, 0, 2); % across columns

SEMOnset   = stdOnset/sqrt(numTrial);  % Standard Error
SEMOffset  = stdOffset/sqrt(numTrial); % Standard Error

% NCIOnset   = 1.96*SEMOnset;  % 95% Confidence Interval
% NCIOffset  = 1.96*SEMOffset; % 95% Confidence Interval

meanAudio = [meanOnset SEMOnset meanOffset SEMOffset];
end

function audioDynamics_Somato = InflationResponse(secTime, secAudioMean)
% audioDynamics_Somato = InflationResponse(secTime, secAudioMean) 
% identifies the relevant pitch contour characteristics that represent how 
% a participant responded to the balloon inflation during vocal production.  
% iR is a structure representing the result variables from studying the 
% inflation response (iR). The prefix letter denotes whether the variable 
% is an index (i), a time (t), or a value (v). 
%
% INPUTS
% secTime:  Vector of time points corresponding to the sectioned data 
%           (numSamp)
% secAudio: 3D mat of sectioned audio (numSamp x numTrial x event)
%           The 1st 3D layer are Onset Sections
%           The 2nd 3D later are Offset Sections
%
% OUTPUTS
% audioDynamics_Somato: Structure of audio dynamics results with the 
%                       following fields:
% -respVarM: Matrix of mean trial values from respVarm (numSamp x 4)

[numSamp, ~] = size(secAudioMean); % Size of the data we are dealing with

ir = initInflationResponseStruct(); % Initialize the structure that handles the variable calculations
ir.time     = secTime;              % Time Interval for the sectioned trials (-0.5->1.0s)
ir.onset    = secAudioMean(:, 1);   % f0 Trace sectioned around pert Onset.

ir.iAtOnset = find(ir.time == 0);
ir.tAtOnset = 0;                     % duh
ir.vAtOnset = ir.onset(ir.iAtOnset); % f0 value at t = 0

ir.iPostOnsetR = find(0 < ir.time & .20 >= ir.time); % Range of indices between t > 0ms and t =< 200ms;
[minOn, minIdx] = min(ir.onset(ir.iPostOnsetR));     % Minimum f0 val within PostOnsetR

% StimMag
ir.iAtMin  = ir.iPostOnsetR(minIdx);       % Indice of the min f0 value following trigger
ir.tAtMin  = ir.time(ir.iAtMin);           % Time at min f0 value following trigger
ir.vAtMin  = minOn;                        % Min f0 value in PostOnsetR
ir.stimMag = abs(ir.vAtMin - ir.vAtOnset); % Distance traveled from onset to min value

% RespMag
ir.tAtRespRange = [0.8 1.0];              % Time Values in the period 800ms to 1000ms after perturbation onset
ir.iAtRespRange = find(ir.time >= ir.tAtRespRange(1) & ir.time<= ir.tAtRespRange(2));                % Last index in section
ir.vAtRespRange = ir.onset(ir.iAtRespRange);    % Time Value when participant 'fully responded' (1.0s)
ir.tAtResp      = mean(ir.tAtRespRange);
ir.vAtResp      = mean(ir.vAtRespRange);
ir.respMag = ir.vAtResp - ir.vAtMin; % Distance traveled from min f0 value to response f0 value

% RespPer
ir.respPer = 100*(ir.respMag/ir.stimMag); % Percent change from stimMag to respMag 

% Add to the audioDynamics struct
respVarM = [ir.tAtMin ir.stimMag ir.respMag ir.respPer];
audioDynamics_Somato.respVarM = respVarM;
% drawInflationResultMetrics(ir, 1, 0, 0); % Generates useful manuscript Fig
end

function ir = initInflationResponseStruct()
% ir = initInflationResponseStruct() initalizes the structure ir which
% organizes the variables tracked for the outcomes measures of an inflation
% response. 
%
% OUTPUTS
% ir: Structure of initalized inflation response variables. The prefix 
% letter denotes whether the variable is an index (i), a time (t), or a 
% value (v)

ir.time     = [];
ir.onset    = [];

ir.iAtOnset = []; % Index where t = 0
ir.tAtOnset = []; % Time at t = 0
ir.vAtOnset = []; % f0 value at t = 0

ir.iPostOnsetR = []; % Range of indices between t = 0ms and t = 200ms;
ir.iAtMin      = []; % Index at min f0 value in PostOnsetR
ir.tAtMin      = []; % Time at min f0 value in PostOnsetR
ir.vAtMin      = []; % Min f0 value in PostOnsetR
ir.stimMag     = []; % ir.vAtMin - ir.vAtOnset ..in a perfect world vAtOnset = 0

ir.iAtResp = []; % Index of f0 value when participant 'fully' responded...right now = last value in section
ir.tAtResp = []; % Time at f0 value when participant 'fully' responded
ir.vAtResp = []; % f0 value when participant 'fully' responded 
ir.respMag = []; % vAtResp - vAtMin   ...distance traveled
ir.respPer = []; % Percent change from stimMag to respMag
end

function audioDynamics_Audio = PitchShiftReflexResponse(secTime, secAudioMean)
% audioDynamics_Audio = PitchShiftReflexResponse(secTime, secAudioMean) 
% identifies the relevant pitch contour characteristics that represent
% how a participant responded to auditory feedback changes during vocal 
% production.
% iR is a structure representing the result variables
% from studying the inflation response (iR). The prefix letter denotes
% whether the variable is an index (i), a time (t), or a value (v). 
%
% INPUTS
% secTime:      Vector of time points corresponding to the sectioned data 
%               (numSamp)
% secAudioMean: 3D mat of sectioned audio (numSamp x numTrial x event)
%               The 1st 3D layer are Onset Sections
%               The 2nd 3D later are Offset Sections
%
% OUTPUTS
% audioDynamics_Audio: Structure of audio dynamics results with the 
%                      following fields:
% -respVarM: Matrix of mean trial values from respVarm (numSamp x 4)

[numSamp, ~] = size(secAudioMean); % Size of the data we are dealing with

ir = initInflationResponseStruct(); % Initialize the structure that handles the variable calculations
ir.time     = secTime;              % Time Interval for the sectioned trials (-0.5->1.0s)
ir.onset    = secAudioMean(:, 1);   % f0 Trace sectioned around pert Onset.

ir.iAtOnset = find(ir.time == 0);
ir.tAtOnset = 0;                     % duh
ir.vAtOnset = ir.onset(ir.iAtOnset); % f0 value at t = 0

ir.iPostOnsetR = find(0 < ir.time & .20 >= ir.time); % Range of indices between t > 0ms and t =< 200ms;
[minOn, minIdx] = min(ir.onset(ir.iPostOnsetR));     % Minimum f0 val within PostOnsetR

% StimMag
ir.iAtMin  = ir.iPostOnsetR(minIdx);       % Indice of the min f0 value following trigger
ir.tAtMin  = ir.time(ir.iAtMin);           % Time at min f0 value following trigger
ir.vAtMin  = minOn;                        % Min f0 value in PostOnsetR
ir.stimMag = abs(-100);                    % Distance traveled from onset to min value (default is 100 cents)

% RespMag
ir.tAtRespRange = [0.8 1.0];              % Time Values in the period 800ms to 1000ms after perturbation onset
ir.iAtRespRange = find(ir.time >= ir.tAtRespRange(1) & ir.time<= ir.tAtRespRange(2));                % Last index in section
ir.vAtRespRange = ir.onset(ir.iAtRespRange);    % Time Value when participant 'fully responded' (1.0s)
ir.tAtResp      = mean(ir.tAtRespRange);
ir.vAtResp      = mean(ir.vAtRespRange);
ir.respMag = ir.vAtResp - ir.vAtMin; % Distance traveled from min f0 value to response f0 value

% RespPer
ir.respPer = 100*(ir.respMag/ir.stimMag); % Percent change from stimMag to respMag 

% Add to the audioDynamics struct
respVarM = [ir.tAtMin ir.stimMag ir.respMag ir.respPer];
audioDynamics_Audio.respVarM = respVarM;
end