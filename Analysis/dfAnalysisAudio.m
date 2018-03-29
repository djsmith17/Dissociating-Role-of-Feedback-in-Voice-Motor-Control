function An = dfAnalysisAudio(dirs, An, AudFlag, varargin)
% An = dfAnalysisAudio(dirs, An, AudFlag, varargin) performs frequency 
% analyses on recorded audio data to measure the pitch contour on a number
% of recorded trials. The actual f0 analysises are done within the
% signalfrequencyanalyses function below, and rest of the script is spent
% converting the f0 traces to cents, and sectioning them around the trigger
% points. Mutli-trial statistics can be perfomed on the set of analyzed
% trials, and inflation response results can be indetified
%
% This function adds to an existing analysis variables structure, which
% means it can use any data set (NIDAQ/Audapter) as long as it has some
% starting variables. 
%
% dirs:    The set of directories we are working in
% An:      Analysis variables used to analyze data set
% AudFlag: Do we want to do analysis on the audio data? Useful for when
%          analyzing just perturbatron data
% iRF:     Inflation response flag. Do we want to analyze the response to
%          inflation?
% f0F:     f0 Flag. Do we want to reanalyze data and receive new pitch
%          traces?
%
% This function expects that the structure An contains the following fields
% -An.curSess
% -An.f0Type
% -An.f0AnaFile
% -An.bTf0b
% -An.sRate
% -An.audioMSvt
% -An.audioHSvt
% -An.expTrigsSvt
% -An.allIdxSvt
% -An.trialTypeSvt

if isempty(varargin)
    iRF    = 0; 
    f0Flag = 0;
elseif length(varargin) == 1
    iRF    = varargin{1};
    f0Flag = 0;
else
    iRF    = varargin{1};
    f0Flag = varargin{2};
end

% Instatiate the variables we intend to create 
An = initAudVar(An);

if AudFlag == 1
    fprintf('\nStarting Pitch Analysis\n')
   
    % Set some frequency analysis variables
    An.numSamp = length(An.audioMSvt);
    fV = setFreqAnalVar(An.sRate, An.numSamp);
    
    % File where to save/find pitch contour analysis
    dirs.audiof0AnalysisFile = fullfile(dirs.SavResultsDir, An.f0AnaFile);
    
    % How do you want to calculate the pitch contour?
    if strcmp(An.f0Type, 'Praat') == 1
        anaFlag = 1;
    else
        anaFlag = 2;
    end
        
    % Pitch contour analysis can be time consuiming. 
    % Flag allows the loading of a previously saved copy for faster reanalysis.
    if exist(dirs.audiof0AnalysisFile, 'file') == 0 || f0Flag == 1

        [f0A.timef0, f0A.audioMf0, f0A.expTrigsR, f0A.etM, f0A.fV] = signalFrequencyAnalysis(dirs, fV, An.audioMSvt, An.expTrigsSvt, An.bTf0b, anaFlag);
        [f0A.timef0, f0A.audioHf0, f0A.expTrigsR, f0A.etH, f0A.fV] = signalFrequencyAnalysis(dirs, fV, An.audioHSvt, An.expTrigsSvt, An.bTf0b, anaFlag);        
        save(dirs.audiof0AnalysisFile, 'f0A')
    else
        load(dirs.audiof0AnalysisFile)
    end

    An.timef0     = f0A.timef0;
    An.expTrigsR  = f0A.expTrigsR;
    An.audioMf0   = f0A.audioMf0;
    An.audioHf0   = f0A.audioHf0;
    An.etMH       = f0A.etM + f0A.etH; % Minutes
    An.fV         = f0A.fV;

    % Smooth the f0 data
    An.audioMf0S   = smoothf0(An.audioMf0);
    An.audioHf0S   = smoothf0(An.audioHf0);

    % Section Audio with all trials...before parsing, and post-processing
    [An.secTime, An.audioMf0SecAll] = sectionData(An.timef0, An.audioMf0S, An.expTrigsR);
    [An.secTime, An.audioHf0SecAll] = sectionData(An.timef0, An.audioHf0S, An.expTrigsR);
   
    % Find the value of f0 during the perPert period for each trial
    prePert       = (An.secTime <= 0); % SecTime is aligned for SecTime = 0 to be Onset of pert
    An.trialf0b   = mean(An.audioMf0SecAll(prePert,:,1),1); % Per-trial baseline f0   
    An.f0b        = mean(An.trialf0b);                      % Mean trial baseline f0
    
    % Normalize f0 traces by individual f0b and convert to cents
    An.audioMf0_norm = normf0(An.audioMf0S, An.trialf0b);
    An.audioHf0_norm = normf0(An.audioHf0S, An.trialf0b);
    
    %Find troublesome trials and remove
    [An.svf0Idx, An.expTrigsf0Sv, An.pertf0Idx, An.contf0Idx] = audioPostProcessing(An);
    
    An.audioMf0sv      = An.audioMf0_norm(:, An.svf0Idx);
    An.audioHf0sv      = An.audioHf0_norm(:, An.svf0Idx); 
    An.numTrialsPP     = length(An.svf0Idx);
    An.numPertTrialsPP = length(An.pertf0Idx);
    An.numContTrialsPP = length(An.contf0Idx);

    %Find the Perturbed Trials
    An.pertTrigsR  = An.expTrigsf0Sv(An.pertf0Idx,:);
    An.audioMf0p = parseTrialTypes(An.audioMf0sv, An.pertf0Idx);
    An.audioHf0p = parseTrialTypes(An.audioHf0sv, An.pertf0Idx);
    An.contTrigsR  = An.expTrigsf0Sv(An.contf0Idx,:);
    An.audioMf0c = parseTrialTypes(An.audioMf0sv, An.contf0Idx);
    An.audioHf0c = parseTrialTypes(An.audioHf0sv, An.contf0Idx);
    
    %Section the data around onset and offset
    [An.secTime, An.audioMf0_Secp] = sectionData(An.timef0, An.audioMf0p, An.pertTrigsR);
    [An.secTime, An.audioHf0_Secp] = sectionData(An.timef0, An.audioHf0p, An.pertTrigsR);
    [An.secTime, An.audioMf0_Secc] = sectionData(An.timef0, An.audioMf0c, An.contTrigsR);
    [An.secTime, An.audioHf0_Secc] = sectionData(An.timef0, An.audioHf0c, An.contTrigsR);

    %Mean around the onset and offset
    An.audioMf0_meanp = meanAudioData(An.audioMf0_Secp);
    An.audioHf0_meanp = meanAudioData(An.audioHf0_Secp);
    An.audioMf0_meanc = meanAudioData(An.audioMf0_Secc);
    An.audioHf0_meanc = meanAudioData(An.audioHf0_Secc); 

    %The Inflation Response
    if iRF == 1
        [An.respVar, An.respVarM, An.respVarSD, An.InflaStimVar] = InflationResponse(An.secTime, An.audioMf0_Secp);
    end
end
end

function An = initAudVar(An)
%Initialize some variables to keep track of them

An.fV             = []; %frequency analysis variables

An.timef0         = []; %time vector of audio samples recorded
An.fsA            = []; %sampling rate of audio samples
An.audioMf0       = []; %Raw Microphone Audio Data
An.audioHf0       = []; %Raw Headphone Audio Data
An.expTrigsR      = [];
An.etMH           = [];
An.audioMf0S      = [];
An.audioHf0S      = [];
An.trialf0b       = []; %Per Trial calculated f0
An.f0b            = []; %Average trial f0
An.audioMf0_norm  = []; %Normalized mic data
An.audioHf0_norm  = []; %Normalized head data

An.svf0Idx        = [];
An.expTrigsf0Sv   = [];
An.pertf0Idx      = [];
An.contf0Idx      = [];
An.audioMf0sv     = [];
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

An.secTime        = [];
An.audioMf0_Secp  = []; 
An.audioHf0_Secp  = [];
An.audioMf0_Secc  = []; 
An.audioHf0_Secc  = [];

An.audioMf0_meanp = []; 
An.audioHf0_meanp = [];
An.audioMf0_meanc = []; 
An.audioHf0_meanc = [];

An.respVar        = []; 
An.respVarM       = [];
An.respVarSD      = [];
An.InflaStimVar   = [];
end

function fV = setFreqAnalVar(sRate, numSamp)

%Identify a few analysis varaibles
fV.sRate      = sRate;
fV.numSamp    = numSamp;
fV.time       = (0:1/fV.sRate:(numSamp-1)/fV.sRate)'; %Time vector for full mic

fV.freqCutOff = 400;
fV.win        = 0.010;      % seconds
fV.fsA        = 1/fV.win;
fV.winP       = fV.win*sRate;
fV.pOV        = 0.60;       % 80% overlap
fV.tStepP     = round(fV.winP*(1-fV.pOV));

fV.trialWin = round(1:fV.tStepP:(fV.numSamp-fV.winP)); % Window start frames based on length of voice onset mic
fV.numWin   = length(fV.trialWin);                     % Number of windows based on WinSt

fV.roundFact = fV.sRate/fV.tStepP;
fV.winHalf   = fV.win/2;
end

function [timef0, audiof0, expTrigsR, elapsed_time, fV] = signalFrequencyAnalysis(dirs, fV, audio, expTrig, bTf0b, flag)
ET = tic;
[~, numTrial] = size(audio);

fs = fV.sRate;

if flag == 1
    fV.win       = 0.005;
    fV.winP      = fV.win*fs;
    fV.pOV       = 0.00;
    fV.tStepP    = round(fV.winP*(1-fV.pOV));
    fV.roundFact = fV.sRate/fV.tStepP;
    fV.winHalf   = 1.0*fV.win;
    
    [timef0, audiof0, fsA] = dfCalcf0Praat(dirs, audio, fs, bTf0b);
else
    audiof0 = [];
    for j = 1:numTrial %Trial by Trial             
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

elapsed_time = toc(ET)/60;
% fprintf('\nElapsed Time: %f (min)\n', elapsed_time)
end

function f0Win = simpleAutoCorr(voice, fV)
%Simple version of an autocorrelation for finding pitch

fs            = fV.sRate;
win           = fV.winP;
[autoC, lags] = autocorr(voice, win-1);
[pks, pkInd]  = findpeaks(autoC);
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
[~, numTrial] = size(audio);

audioS = [];
for ii = 1:numTrial
    audioSmooth = smooth(audio(:,ii), 10);   % 10 sample length smoothing
    audioS      = cat(2, audioS, audioSmooth);
end
end

function audio_norm = normf0(audio, f0b)
% audio_norm = normf0(audio, f0b) takes a matrix of audio signals (audio) 
% of size numSamp x numTrial and normalizes each trial by the f0 caluclated 
% for that trial which are stored in the vector f0b (numTrial x 1)

[~, numTrial] = size(audio);

audio_norm = [];
for ii = 1:numTrial
    audio_trial = 1200*log2(audio(:,ii)./f0b(ii));
    audio_norm  = cat(2, audio_norm, audio_trial);
end
end

function signalParse = parseTrialTypes(signal, idx)
% signalParse = parseTrialTypes(signal, idx) parses individual trials out 
% of a large matrix of recordings of size numSamp x numTrial. 
% (idx) is a vector of the indices to parse out.
% Why did you make this a function? Get over it.

signalParse = signal(:, idx);
end

function [svf0Idx, expTrigsf0Sv, pertf0Idx, contf0Idx] = audioPostProcessing(An)
% [svf0Idx, expTrigsf0Sv, pertf0Idx, contf0Idx] = audioPostProcessing(An)
% checks for odd (innaproriately large or small) normalized f0 values in
% whole f0 traces as a result of the frequency analyses. 
% This throws away trials any such odd trials and prints a statement about
% which trials were thrown out.
%
% This function returns
% svf0Idx: The indices of the saved trials compared against the full set
% expTrigsf0Sv: The trigs of the saved trials
% pertf0Idx: The indices of svf0Idx that are perturbed trials
% contf0Idx: The indices of svf0Idx that are control trials

curSess      = An.curSess;
svIdx        = An.allIdxSvt;
trialTypeSvt = An.trialTypeSvt;

time        = An.timef0;
audioNormM  = An.audioMf0_norm;
expTrigs    = An.expTrigsR;

[~, numTrialType] = size(audioNormM);

timeInd = (time > 0.5 & time < 3.5);

svf0Idx      = [];
expTrigsf0Sv = [];
contf0Idx    = [];
pertf0Idx    = [];
svF = 0;
for ii = 1:numTrialType
    
    mic     = audioNormM(timeInd, ii);
    expTrig = expTrigs(ii, :);
    svIdc   = svIdx(ii);
    
    ind = find(mic >= 500 | mic <=  -500);
    
    if ~isempty(ind)
        if trialTypeSvt(ii) == 0
            type = 'Cont';
        else
            type = 'Pert';      
        end       
        fprintf('Threw away %s Trial %d (%s)\n', curSess, svIdc, type)
        
    else
        svF = svF + 1;
        
        svf0Idx      = cat(1, svf0Idx, ii);
        expTrigsf0Sv = cat(1, expTrigsf0Sv, expTrig);
        if trialTypeSvt(ii) == 0
            contf0Idx  = cat(1, contf0Idx, svF);
        else
            pertf0Idx  = cat(1, pertf0Idx, svF);       
        end        
    end
end
end

function [secTime, secSigs] = sectionData(time, sigs, trigs)
% [secTime, secSigs] = sectionData(time, sigs, trigs) sections
% time series data around important points in time.
% 
% time:  Vector of time points (numSamp)
% sigs:  Matrix of signals to be sectioned (numSamp x numTrial)
% trigs: Onset and Offset time tiggers (numTrial x 2)
%
% secTime: Vector of time points corresponding to the sectioned window (numSampSec)
% secSigs: 3D mat of sectioned sigs (numSampSec x numTrial x event)
%          The 1st 3D layer are Onset Sections
%          The 2nd 3D later are Offset Sections

[~, numTrial] = size(sigs);
preEve  = 0.5; posEve = 1.0;

secSigs    = [];
OnsetSecs  = [];
OffsetSecs = [];
for ii = 1:numTrial
    OnsetT   = trigs(ii, 1); % Onset time
    OffsetT  = trigs(ii, 2); % Offset time
    
    OnsetTSt = round(OnsetT - preEve, 3);   % PreOnset time, rounded to nearest ms
    OnsetTSp = round(OnsetT + posEve, 3);   % PostOnset time, rounded to nearest ms
    OnsetSpan = time >= OnsetTSt & time <= OnsetTSp; % Indices corresponding to Onset period
    
    OffsetTSt = round(OffsetT - preEve, 3); % PreOffset time, rounded to nearest ms
    OffsetTSp = round(OffsetT + posEve, 3); % PostOffset time, rounded to nearest ms
    OffsetSpan = time >= OffsetTSt & time <= OffsetTSp; % Indices corresponding to Offset period
        
    OnsetSec  = sigs(OnsetSpan, ii);  % Data sectioned around Onset
    OffsetSec = sigs(OffsetSpan, ii); % Data sectioned around Offset
    
    OnsetSecs  = cat(2, OnsetSecs, OnsetSec);   % Sectioned signal onsets concatenated
    OffsetSecs = cat(2, OffsetSecs, OffsetSec); % Sectioned signal offsets concatenated
end

numSampSec = length(OnsetSec); % number of samples in sectioned signals

secTime = linspace(-preEve, posEve, numSampSec); % time vector correspnding to the sectioned signals
secSigs(:,:,1) = OnsetSecs;  % 1st 3D layer
secSigs(:,:,2) = OffsetSecs; % 2nd 3D layer
end

function y = round2matchfs(x, rFact, winHalf)
% y = round2matchfs(x, rFact, winHalf) rounds the value of x to a value y,
% which is the closet number to x that matches the time step for a given
% analysis windowing methods. winHalf is the length of the window that
% represents the amount of window, the time point corresponds to. rFact is
% the sampling rate divided by the step size.
%
% x can be a single number, a vector of numbers or a matrix of numbers.

y = round((x-winHalf).*rFact)./rFact + winHalf;
end

function meanAudio = meanAudioData(secAudio)
% Some simple statistics on the sectioned audio data. 
% meanAudio is a vector containing the following information
% meanAudio(1) = mean Onset pitch contour
% meanAudio(2) = 95% CI of the mean Onset Pitch Contour
% meanAudio(3) = mean Offset pitch contour
% meanAudio(4) = 95% CI of the mean Offset Pitch Contour

OnsetSecs  = secAudio(:,:,1);
OffsetSecs = secAudio(:,:,2);
[~, numTrial] = size(OnsetSecs);

meanOnset  = mean(OnsetSecs, 2);  % across columns
meanOffset = mean(OffsetSecs, 2); % across columns

stdOnset   = std(OnsetSecs, 0, 2);  % across columns
stdOffset  = std(OffsetSecs, 0, 2); % across columns

SEMOnset   = stdOnset/sqrt(numTrial);  % Standard Error
SEMOffset  = stdOffset/sqrt(numTrial); % Standard Error

NCIOnset   = 1.96*SEMOnset;  % 95% Confidence Interval
NCIOffset  = 1.96*SEMOffset; % 95% Confidence Interval

meanAudio = [meanOnset NCIOnset meanOffset NCIOffset];
end

function [respVar, respVarM, respVarSD, InflaStimVar] = InflationResponse(secTime, secAudio)
% [respVar, respVarm, respVarSD, InflaStimVar] = InflationResponse(secTime, secAudio)
% Identifies the relevant pitch contour characteristics that are important
% for deciding how a participant responded to the inflation of the balloon
% during production. iR is a structure representing the result variables
% from studying the inflation response (iR). 
%
% secTime:  Vector of time points corresponding to the sectioned data (numSamp)
% secAudio: 3D mat of sectioned audio (numSamp x numTrial x event)
%           The 1st 3D layer are Onset Sections
%           The 2nd 3D later are Offset Sections
%
% respVar: Matrix of per trial iR results (numTrial x 4)
%          respVar(:,1) = Time of the minimum f0 value in the sec
%          respVar(:,2) = Minimum f0 value in sec (stim magnitude)
%          respVar(:,3) = Value of f0 at end of sec (response magnitude)
%          respVar(:,4) = ABS percent change of stim and response
%          magnitude (response percentage)
% respVarM:    Vector of mean trial values from respVarm (1x4)
% respVarSD:   Vector of standard deviation of the trial values from respVar (1x4)
% InflaSimVar: Values of the mean time at stim magnitude and the mean stim magnitude

[numSamp, numTrial, ~] = size(secAudio); % Size of the data we are dealing with

ir.numSamp  = numSamp;
ir.numTrial = numTrial;
ir.time     = secTime;

ir.iAtOnset = find(secTime == 0); % Index where t = 0
ir.tAtOnset = 0;                  % Time at t = 0 ...duh
ir.vAtOnset = [];                 % f0 value at t = 0

ir.iPostOnsetR = find(0 <= secTime & .20 >= secTime); % Range of indices between t = 0ms and t = 200ms;
ir.iAtMin  = [];                  % Index at min f0 value in PostOnsetR
ir.tAtMin  = [];                  % Time at min f0 value in PostOnsetR
ir.vAtMin  = [];                  % Min f0 value in PostOnsetR
ir.stimMag = [];                  % ir.vAtMin - ir.vAtOnset ..in a perfect world vAtOnset = 0

ir.iAtResp = ir.numSamp;          % Index of f0 value when participant 'fully' responded...right now = last value in section
ir.tAtResp = ir.time(ir.numSamp); % Time at f0 value when participant 'fully' responded
ir.vAtResp = [];                  % f0 value when participant 'fully' responded 
ir.respMag = [];                  % vAtResp - vAtMin   ...distance traveled
ir.respPer = [];                  % Percent change from stimMag to respMag

% Variables to be concatenated and saved as outputs 
tAtMin  = []; stimMag = [];
respMag = []; respPer = [];
for i = 1:numTrial
    onset = secAudio(:,i,1); % Go trial by trial; First 3D layer is Onset
    ir.vAtOnset = onset(ir.iAtOnset); % f0 value at t = 0

    [minOn, minIdx] = min(onset(ir.iPostOnsetR)); % Minimum f0 in PostOnsetR
    ir.iAtMin = ir.iPostOnsetR(minIdx);           % Indice of the min f0 value
    ir.tAtMin = ir.time(ir.iAtMin);               % Time at min f0 value in PostOnsetR
    ir.vAtMin = minOn;                            % Min f0 value in PostOnsetR
    ir.stimMag = ir.vAtMin - ir.vAtOnset;         % Distance traveled from onset to min value
    
    ir.vAtResp = onset(ir.iAtResp);               % f0 value when participant 'fully' responded 
    ir.respMag = ir.vAtResp - ir.vAtMin;          % Distance traveled from min f0 value to response f0 value
    
    ir.respPer = 100*(ir.respMag/abs(ir.stimMag));% Percent change from stimMag to respMag 
    
    if ir.stimMag == 0
        ir.respPer = 0.0;
    end
    
    % Concatenate the results from this trial 
    tAtMin   = cat(1, tAtMin, ir.tAtMin);
    stimMag  = cat(1, stimMag, ir.stimMag); 
    respMag  = cat(1, respMag, ir.respMag); 
    respPer  = cat(1, respPer, ir.respPer);
end

% Organize the results 
respVar   = [tAtMin stimMag respMag respPer];
respVarM  = mean(respVar, 1);
respVarSD = std(respVar, 0, 1);

InflaStimVar = [respVarM(1) respVarM(2)];
end