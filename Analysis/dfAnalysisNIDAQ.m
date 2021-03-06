function [niAn, niRes] = dfAnalysisNIDAQ(dirs, expParam, DAQin, AudFlag, aDF, PresFlag)
% [niAn, niRes] = dfAnalysisNIDAQ(dirs, expParam, DAQin, AudFlag, aDF, PresFlag)
% This function organizes and analyzes the raw NIDAQ recordings taken from
% sensors and audio devices during experiments measuring changes in f0. 
% Specifically this function is responsible for cataloging the
% timing dynamics of the perturbation trigger signals, as well as measuring
% the sensor dynamics of the sensors that verify appropriate stimulus
% levels
% 
% dirs:     The set of directories we are currently working in 
% expParam: The experimental parameters of the recorded experiment
% DAQin:    Raw NIDAQ data structure
% AudFlag:  Flag to check if analyses of audio data should be performed
% aDF:      Audio Dynamics Flag. Analyze changes in f0 following triggers
% PresFlag: Flag to check if analyses of pressure data should be performed
%
% niAn:  Structure of all variables used to analyze NIDAQ data
% niRes: Structure of necessary vars for pooled-analyses, stats, & figures
%
% This function calls the following functions
% -dfAnalysisAudio.m
%
% See below for the following sub-functions:
% -initNIDAQAnalysisStruct
% -initSensorDynamicsStruct
% -convertPressureSensor
% -analyzeSensorDynamics
%
% Requires the Signal Processing Toolbox, Image Processing Toolbox

%Identify some starting variables
niAn = initNIDAQAnalysisStruct();

niAn.expType   = expParam.expType;
niAn.subject   = expParam.subject;
niAn.run       = expParam.run;
niAn.curSess   = expParam.curSess;

niAn.f0AnaFile = [niAn.subject niAn.run 'f0Analysis.mat'];
niAn.gender    = expParam.gender;
niAn.AudFB     = expParam.AudFB;
niAn.AudFBSw   = expParam.AudFBSw;

if isfield(expParam, 'f0b')
    niAn.bTf0b     = expParam.f0b;
else
    niAn.bTf0b     = 100;
end

if isfield(expParam, 'balloon')
    niAn.balloon = expParam.balloon;
else
    niAn.balloon = 'N/A';
end

if isfield(expParam, 'sensorPType')
    niAn.sensorPType = expParam.sensorPType;
else
    niAn.sensorPType = 'Five';
end

fprintf('Starting NIDAQ Analysis for %s, %s with f0 of %0.2f Hz\n', niAn.subject, niAn.run, niAn.bTf0b)

[r, c, n]      = size(DAQin);
niAn.sRate     = expParam.sRateQ;       % Sampling Rate of the NIDAQ
niAn.numCh     = c;                     % Number of Channels recorded
niAn.numSamp   = r;                     % Number of Samples recorded
niAn.numTrial  = n;                     % Number of Trials recorded
niAn.trialLen  = niAn.numSamp/niAn.sRate;
niAn.trialType = expParam.trialType;    % Control (0), Perturbed (1)
niAn.expTrigs  = expParam.trigs(:,:,1); % Trigger Onset and Offset (Time) (all recorded trials)

% Find all the perturbed trials
[niAn.ContTrials, niAn.contIdx] = find(niAn.trialType == 0);
[niAn.PertTrials, niAn.pertIdx] = find(niAn.trialType == 1);
niAn.numContTrials = sum(niAn.ContTrials);
niAn.numPertTrials = sum(niAn.PertTrials);

%Unpack the NIDAQ raw data set, a 3D matrix (numSamp, numCh, numTrial)
niAn.time     = linspace(0, niAn.trialLen, niAn.numSamp)';
niAn.pertSig  = squeeze(DAQin(:,1,:)); % Perturbatron Signal (Pert)
niAn.sensorFC = squeeze(DAQin(:,2,:)); % Force Sensor Collar (FC)
niAn.sensorFN = squeeze(DAQin(:,3,:)); % Force Sensor Neck   (FN)
niAn.sensorP  = squeeze(DAQin(:,4,:)); % Pressure            (P)
niAn.audioM   = squeeze(DAQin(:,5,:)); % Microphone Signal   (M)
niAn.audioH   = squeeze(DAQin(:,6,:)); % Headphone Signal    (H)
niAn.sensorO  = squeeze(DAQin(:,7,:)); % Optical Trigger Box (O)

% Convert the measured Pressure Sensor Voltage (V) to Pressure (psi)
niAn.sensorFNz = convertPressureSensor(niAn.sensorFN, niAn.sensorPType);
niAn.sensorPz  = convertPressureSensor(niAn.sensorP, niAn.sensorPType);

% Parse out the perturbed trials
niAn.pertSig_p  = parseTrialTypes(niAn.pertSig, niAn.pertIdx);   % Only Perturbed Trials
niAn.sensorP_p  = parseTrialTypes(niAn.sensorPz, niAn.pertIdx);  % Only Perturbed Trials
niAn.sensorFC_p = parseTrialTypes(niAn.sensorFC, niAn.pertIdx);  % Only Perturbed Trials
niAn.sensorFN_p = parseTrialTypes(niAn.sensorFNz, niAn.pertIdx); % Only Perturbed Trials

%Find Rising and Falling Edges of sensor signals: Onset and Offset TRIGGERS
[niAn.pertTrig, niAn.idxPert] = findPertTrigs(niAn.time, niAn.pertSig_p);

% What was the delay between the intended (code) onset/offset triggers
% and actual (measured) onset/offset triggers
OnsetTriggerLags  = round((niAn.pertTrig(:, 1) - niAn.expTrigs(niAn.pertIdx, 1)), 3);
OffsetTriggerLags = round((niAn.pertTrig(:, 2) - niAn.expTrigs(niAn.pertIdx, 2)), 3);
niAn.TriggerLag   = [OnsetTriggerLags, OffsetTriggerLags];

niAn.alignResponseTriggers = niAn.expTrigs;

niAn.pertSD.TrigTime = niAn.pertTrig;
niAn.pertSD.TrigIdx  = niAn.idxPert;

niAn.presSD.time   = niAn.time;
niAn.presSD.sensor = niAn.sensorP_p;
niAn.presSD.fs     = niAn.sRate;

niAn.fSNSD.time   = niAn.time;
niAn.fSNSD.sensor = niAn.sensorFN_p;
niAn.fSNSD.fs     = niAn.sRate;

if PresFlag == 1 && niAn.numPertTrials > 0
    % Set PresFlag = 1 if pressure dynamics are worth looking investigating
    % Observe dynamics of the pressure sensor and save pert-onset aligned
    % recordings. Also saving a set of data set that is sectioned around 
    % the onset and offset of the perturbation period.
    
    % Pressure Sensor Dynamics
    [niAn.presSD] = analyzeSensorDynamics(niAn.presSD, niAn.pertSD);
    
%     SD = StepFunctionDynamics(niAn.time, niAn.sensorP_p, niAn.sRate, niAn.idxPert, niAn.pertTrig, niAn.curSess);
%     SD.drawAllTrial;
    
    [niAn.fSNSD] = analyzeSensorDynamics(niAn.fSNSD, niAn.pertSD);
end

%The Audio Analysis
niAn = dfAnalysisAudio(dirs, niAn, AudFlag, aDF);
    
lims  = identifyLimits(niAn);
niRes = packResults(niAn, lims);
end

function niAn = initNIDAQAnalysisStruct()

niAn.AnaType   = 'NIDAQ';
niAn.expType   = [];
niAn.subject   = [];
niAn.run       = [];
niAn.curSess   = [];
niAn.f0Type    = 'Praat';
niAn.f0AnaFile = [];
niAn.gender    = [];
niAn.AudFB     = [];
niAn.AudFBSw   = [];
niAn.bTf0b     = [];

niAn.balloon     = [];
niAn.sensorPType = [];

niAn.sRate     = [];
niAn.numCh     = [];
niAn.numSamp   = [];
niAn.numTrial  = [];
niAn.trialLen  = [];
niAn.trialType = [];
niAn.expTrigs  = [];

niAn.ContTrials    = [];
niAn.contIdx       = [];
niAn.PertTrials    = [];
niAn.pertIdx       = [];
niAn.numContTrials = [];
niAn.numPertTrials = [];

niAn.time     = [];
niAn.pertSig  = [];
niAn.sensorFC = [];
niAn.sensorFN = [];
niAn.sensorP  = [];
niAn.audioM   = [];
niAn.audioH   = [];
niAn.sensorO  = [];

niAn.sensorPz  = [];
niAn.sensorFCz = [];
niAn.sensorFNz = [];

niAn.pertSig_p  = [];
niAn.sensorP_p  = [];
niAn.sensorFC_p = [];
niAn.sensorFN_p = [];

niAn.pertTrig   = [];
niAn.idxPert    = [];
niAn.TriggerLag = [];
niAn.alignResponseTriggers = [];

niAn.pertSD = initSensorDynamicsStruct();
niAn.presSD = initSensorDynamicsStruct();
niAn.fSCSD  = initSensorDynamicsStruct();
niAn.fSNSD  = initSensorDynamicsStruct();
end

function SD = initSensorDynamicsStruct()

SD.time   = [];
SD.sensor = [];
SD.fs     = [];

SD.pertIdx  = [];
SD.pertTime = [];

% Per Trial sensor index and time value of rising edge start and falling
% edge end (expecting ~step function)
SD.TrigIdx  = [];
SD.TrigTime = [];

% Per Trial lag time of trigger to start of rising edge
SD.lagTimes  = [];
SD.lagTimeM  = [];
SD.lagTimeSE = [];

% Per Trial rise time of rising edge (step function)
SD.riseTimes  = [];
SD.riseTimeM  = [];
SD.riseTimeSE = [];

% Per trial value at Onset/Offset (pressure, voltage, etc)
SD.OnOffVal   = [];
SD.OnOffValM  = [];
SD.OnOffValSE = [];

% Per trial loss in value (pressure, voltage, etc)
SD.pTrialLoss   = [];
SD.pTrialLossM  = [];
SD.pTrialLossSE = [];

% Sensor recording aligned at perturbation trigger
SD.timeAl   = [];
SD.sensorAl = [];

% Sensor recording sectioned around onset/offset
SD.timeSec    = [];
SD.sensorSec  = [];
SD.sensorSecM = [];
end

function sensorPres = convertPressureSensor(sensorV, sensorType)
% sensorPres = convertPressureSensor(sensorV, sensorType) converts the 
% the recorded voltage from the pressure sensor to the actual pressure
% of the system. The conversion is doing using the transfer function of the
% sensor circuit itself. The sensor being used will have a different Max
% Voltage and therefore a slightly different function. See the attached 
% links for information about the sensors used and their transfer functions

switch sensorType
    case 'Five'
        PMax      = 5;
        PMin      = 0;
        VMax      = 4.5;
        VMin      = 0.5;
        Vsupply   = 5.2;
    case 'Seven'
        PMax      = 7.25;
        PMin      = 0;
        VMax      = 4.5;
        VMin      = 0.5;
        Vsupply   = 5.2;
end

% sensorPres = (sensorV - 0.5)*PMax/4;
sensorPres = PMin + (sensorV - 0.1*Vsupply)*(PMax - PMin)/(0.8*Vsupply);
% 
% m = (PMax - PMin) / (VMax - VMin);
% b = PMin - m*VMin;
% 
% sensorPres = sensorV*m + b; % Convert from voltage to pressure
end

function signalParse = parseTrialTypes(signal, idx)
% signalParse = parseTrialTypes(signal, idx) parses individual trials out 
% of a large matrix of recordings of size numSamp x numTrial. 
% (idx) is a vector of the indices to parse out.
% Why did you make this a function? Get over it. 

signalParse = signal(:, idx);
end

function [trigs, idx] = findPertTrigs(time, sensor)
%findPertTrigs(time, sensor) finds rising and falling edges in sensor
%data. It is expected that these signals will be mostly step functions
[~, numTrial] = size(sensor);

trigs = [];
idx   = [];
for i = 1:numTrial
    
    fDiff = [0; diff(sensor(:,i))];
    fDiff = smooth(fDiff);
    
    threshUp = 0.5*max(fDiff);
    threshDn = 0.5*min(fDiff);
    ups = find(fDiff > threshUp);
    dns = find(fDiff < threshDn);
    
    idxSt = ups(1); 
    idxSp = dns(1);       
    trigSt = round(time(idxSt), 3);
    trigSp = round(time(idxSp), 3);
    
    trigs = cat(1, trigs, [trigSt trigSp]);
    idx   = cat(1, idx, [idxSt idxSp]);
end
end

function [SD] = analyzeSensorDynamics(SD, pertSD)
% Analyzing the dynamics of the sensor during onset and offset of the
% stimulus signal.

time     = SD.time;
sensor   = SD.sensor;
fs       = SD.fs;

SD.pertIdx  = pertSD.TrigIdx;
SD.pertTime = pertSD.TrigTime;

fiveMs = round(0.005*fs); %5ms in points

[numSamp, numTrial] = size(sensor);

for ii = 1:numTrial
    trial = sensor(:,ii);
    
    trialSmoothed = smooth(trial, 100);

    % 1st Derivative (1D) of the recording
    fDiff = 10*[0; diff(trialSmoothed)];

    % Shave off the 5ms at beginning and end
    fDiff(1:fiveMs) = 0;
    fDiff((numSamp-fiveMs):end) = 0;
    
    % Thresholds for detecting edges of (expected) step function
    threshUp = 0.3*max(fDiff);
    threshDn = 0.3*min(fDiff);
    
    %%%Find the Rising Edge%%%
    ups = find(fDiff > threshUp); % Parts of 1D that could include Increasing edge
    StRiseIdx = ups(1);           % Assume first idx of increasing edges is (St)art of rise

    RisingEdgeRange = StRiseIdx:(StRiseIdx + 0.3*fs); % Range Following the StRiseIdx
    [~, idxAtMax] = max(trial(RisingEdgeRange));      % Max value of recording in that range (Assume low->high)
    SpRiseIdx = StRiseIdx + idxAtMax-1;               % Idx of Rise (S)to(p)

    %%%Find the Falling Edge%%%
    followRiseEdge = fDiff(StRiseIdx:end);
    dns = find(followRiseEdge < threshDn); % Parts of 1D that could include decreasing edge                 
    StFallIdx = dns(1)+ StRiseIdx;           % Assume first idx of decreasing edges is (St)art of fall

    FallingEdgeRange = StFallIdx:(StFallIdx + 0.3*fs); % Range Following the StFallIdx
    [~, idxAtMin] = min(trial(FallingEdgeRange));      % Min value of recording in that range (Assume high->low)
    SpFallIdx = StFallIdx + idxAtMin-1;                % Idx of Fall (S)to(p)

    % Convert Indices to Times
    StRiseTime = round(time(StRiseIdx), 3);
    SpRiseTime = round(time(SpRiseIdx), 3);
    StFallTime = round(time(StFallIdx), 3);
    SpFallTime = round(time(SpFallIdx), 3);
    
    lagTimeRise = StRiseTime - SD.pertTime(ii, 1); 
    lagTimeFall = StFallTime - SD.pertTime(ii, 2);
    
    riseTime = SpRiseTime - StRiseTime;
    fallTime = SpFallTime - StFallTime;
    
    % What was the val of the sensor at the end of Rising Edge (Start Plateau)
    RiseVal = trial(SpRiseIdx);
    
    % What was the val of the sensor at the start of the Falling Edge (End Plateau)
    FallVal = trial(StFallIdx);

    SD.TrigIdx  = cat(1, SD.TrigIdx, [StRiseIdx StFallIdx]);
    SD.TrigTime = cat(1, SD.TrigTime, [StRiseTime StFallTime]);
    
    SD.lagTimes  = cat(1, SD.lagTimes, [lagTimeRise lagTimeFall]);
    SD.riseTimes = cat(1, SD.riseTimes, [riseTime fallTime]);
    
    SD.OnOffVal = cat(1, SD.OnOffVal, [RiseVal, FallVal]);
    
%     figure
%     plot(time, trial)
%     hold on
%     plot([StRiseTime StRiseTime], [-5 10], 'g')
%     hold on
%     plot([SpRiseTime SpRiseTime], [-5 10], 'g--')
%     hold on
%     plot([StFallTime StFallTime], [-5 10], 'r')
%     hold on
%     plot([SpFallTime SpFallTime], [-5 10], 'r--')
%     axis([0 4 -0.5 5])
end

%Difference in val at between subsequent trials at SpRiseIdx
SD.pTrialLoss = diff(SD.OnOffVal(:,1));

%Mean and SE lag time
SD.lagTimeM  = mean(SD.lagTimes);
SD.lagTimeSE = (std(SD.lagTimes))/sqrt(numTrial);

%Mean and SE rise time
SD.riseTimeM  = mean(SD.riseTimes);
SD.riseTimeSE = (std(SD.riseTimes))/sqrt(numTrial);

%Mean and SE Val at SpRiseIdx
SD.OnOffValM  = mean(SD.OnOffVal);
SD.OnOffValSE = (std(SD.OnOffVal))/sqrt(numTrial);

%Mean and SE trial-trial difference of val at SpRiseIdx
SD.pTrialLossM  = mean(SD.pTrialLoss);
SD.pTrialLossSE = (std(SD.pTrialLoss))/sqrt(numTrial-1);

% Round and Covert to standard units
SD.lagTimeM  = round(1000*SD.lagTimeM, 3); % Round and Convert s -> ms
SD.lagTimeSE = round(1000*SD.lagTimeSE, 3); % Round and Convert s -> ms

SD.riseTimeM  = round(1000*SD.riseTimeM, 3); % Round and Convert s -> ms
SD.riseTimeSE = round(1000*SD.riseTimeSE, 3); % Round and Convert s -> ms

SD.OnOffValM  = round(SD.OnOffValM, 3); % Round
SD.OnOffValSE = round(SD.OnOffValSE, 3); % Round

SD.pTrialLossM  = round(SD.pTrialLossM, 3); % Round
SD.pTrialLossSE = round(SD.pTrialLossSE, 3); % Round

%%%%%%%%%%%%%
[SD.timeAl, SD.sensorAl] = alignSensorData(sensor, fs, SD.TrigIdx);

[SD.timeSec, SD.sensorSec] = sectionData(sensor, fs, SD.TrigIdx);
SD.sensorSecM              = meanSensorData(SD.sensorSec);   
end

function [timeAl, sensorAl] = alignSensorData(sensor, fs, idx)
% [timeAl, sensorAl]  = alignSensorData(sensor, fs, idx) sections 
% multi-trial sensor data about individual trial trigger points. 
% Each sectioned trial is of equal length, and includes equal lengths of 
% data on either side of the trigger point. 
% The sectioned trials are then concatenated into a matrix, which are
% aligned the trigger point of each trial. 
%
% sensor: recorded sensor data (numSamp x numTrial)
% fs:     sampling rate of sensor data
% idx:    Onset and Offset trigger POINTS (numTrial, 2)
%
% timeAl:   vector of time points corresponding to the sectioned data
% sensorAl: sectioned and aligned sensor data 

[~, numTrial] = size(sensor);
preEve = 0.5; % time preEvent Seconds 
posEve = 2.0; % time posEvent Seconds

% At the moment only aligning by Onset. 
% This could eventually become an input to toggle between Onset/Offset
OnsetTrigs = idx(:, 1); 

sensorAl = [];
for ii = 1:numTrial
    St = OnsetTrigs(ii) - fs*preEve;  % Points preEvent (trigger)
    Sp = OnsetTrigs(ii) + fs*posEve;  % Points posEvent (trigger)
    
    sensorSec = sensor(St:Sp, ii); % Grab St:Sp around the trigger point for this trial
    
    sensorAl = cat(2, sensorAl, sensorSec);
end
per     = 1/fs;

% Time Vector of the sectioned data
timeAl = (-preEve:per:posEve)';
end

function [secTime, secSigs] = sectionData(sigs, fs, trigs)
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
preEveT = 0.5; posEveT = 1.0;
preEve  = preEveT*fs; posEve = posEveT*fs;

secSigs    = [];
OnsetSecs  = [];
OffsetSecs = [];
if numTrial > 0
    for ii = 1:numTrial
        thisSig = sigs(:, ii);
        thisSig = padarray(thisSig,1000,0,'post');
        
        Onset   = trigs(ii, 1); % Onset point
        Offset  = trigs(ii, 2); % Offset point

        OnsetPre  = Onset - preEve;   % PreOnset point
        OnsetPos  = Onset + posEve;   % PostOnset point
        OnsetSpan = OnsetPre:OnsetPos; % Indices corresponding to Onset period

        OffsetPre  = Offset - preEve;   % PreOnset point
        OffsetPos  = Offset + posEve;   % PostOnset point
        OffsetSpan = OffsetPre:OffsetPos; % Indices corresponding to Onset period

        OnsetSec  = thisSig(OnsetSpan);  % Data sectioned around Onset
        OffsetSec = thisSig(OffsetSpan); % Data sectioned around Offset

        OnsetSecs  = cat(2, OnsetSecs, OnsetSec);   % Sectioned signal onsets concatenated
        OffsetSecs = cat(2, OffsetSecs, OffsetSec); % Sectioned signal offsets concatenated
    end
    [numSampSec, ~] = size(OnsetSecs); % number of samples in sectioned signals
else
    numSampSec = 12000;
end

secTime = linspace(-preEveT, posEveT, numSampSec); % time vector correspnding to the sectioned signals
secSigs(:,:,1) = OnsetSecs;  % 1st 3D layer
secSigs(:,:,2) = OffsetSecs; % 2nd 3D layer
end

function meanAudio = meanSensorData(secAudio)
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

% NCIOnset   = 1.96*SEMOnset;  % 95% Confidence Interval
% NCIOffset  = 1.96*SEMOffset; % 95% Confidence Interval

meanAudio = [meanOnset SEMOnset meanOffset SEMOffset];
end

function lims = identifyLimits(An)
% identifyLimits(An) calculates limits of analyzed data so that the limits
% are dynamic and fit the data being shown. 
%
% lims is a structure of the resultant limits to be used in plotting. 
% This function is redundant between a few different functions, and might
% eventually become its own function

%Full Inidividual Trials: Pressure Sensor
maxPres = max(max(An.sensorP_p));
minPres = min(min(An.sensorP_p));
uLP = maxPres + 0.2;
lLP = minPres - 0.2;

lims.pressure   = [0 4 lLP uLP];

%Aligned Pressure Data
lims.pressureAl = [-1 2.5 lLP uLP];

%Mean Pressure Data
lims.pressureMean = [-0.5 1.0 lLP uLP];

%Full Individual Trials: Force Sensors
lims.force      = [0 4 1 5];

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

function res = packResults(niAn, lims)
% packResults(niAn, lims) takes the results of the analysis and packages
% the important variables into a new structure that have common names
% between other analysis methods. This makes it easy to switch back and
% forth between different result structures for plotting. 

% Identifying Subject Information
res.subject = niAn.subject;
res.run     = niAn.run;
res.curSess = niAn.curSess;

% Identifying Experimental Settings
res.expType = niAn.expType;
res.AudFB   = niAn.AudFB;
res.balloon = niAn.balloon;

res.numTrials     = niAn.numTrial;
res.numContTrials = niAn.numContTrials;
res.numPertTrials = niAn.numPertTrials;
res.contIdx       = niAn.contIdx;
res.pertIdx       = niAn.pertIdx;
res.pertTrig      = niAn.pertTrig;

% How would you like the audio to be aligned? Against the Trigger Onset?
% Against the pressure onset?
res.alignResponseTriggers = niAn.alignResponseTriggers;

res.timeS         = niAn.time;
% Pressure Results
res.sensorP       = niAn.sensorP_p; % Individual Processed Perturbed trials. 
res.presSD        = niAn.presSD;    % Sensor Dynamics Structure
res.limitsP       = lims.pressure;  % Limits for collection of individual trials

% Force Sensor (Neck) Results
res.fSNSD         = niAn.fSNSD;     % Sensor Dynamics Structure

% Limits for Aligned and Meaned Pressure Recordings
res.limitsPAl   = lims.pressureAl;
res.limitsPMean = lims.pressureMean;

% Audio f0 analysis
res.timef0          = niAn.timef0;
res.f0b             = niAn.trialf0M;

res.numContTrialsPP = niAn.numContTrialsPP;
res.numPertTrialsPP = niAn.numPertTrialsPP;
res.pertTrigPP      = niAn.pertTrigsR;

%Full Individual Trials: Mic/Head f0 Trace 
res.audioMf0TrialPert = niAn.audioMf0p;
res.audioMf0TrialCont = niAn.audioMf0c;
res.audioHf0TrialPert = niAn.audioHf0p;
res.audioHf0TrialCont = niAn.audioHf0c;
res.limitsA           = lims.audioM;
res.limitsAudRes      = lims.audioAudRespMH;

%Sections Trials: Mic/Head f0
res.secTime          = niAn.secTime;
res.audioMf0SecPert  = niAn.audioMf0_Secp;
res.audioMf0SecCont  = niAn.audioMf0_Secc;
res.audioHf0SecPert  = niAn.audioHf0_Secp;
res.audioHf0SecCont  = niAn.audioHf0_Secc;

%Mean Sectioned Trials: Mic/Head f0 Trace 
res.audioMf0MeanPert = niAn.audioMf0_meanp; % [MeanSigOn 90%CI MeanSigOff 90%CI]
res.audioMf0MeanCont = niAn.audioMf0_meanc;
res.audioHf0MeanPert = niAn.audioHf0_meanp;
res.audioHf0MeanCont = niAn.audioHf0_meanc;
res.limitsAmean      = lims.audioMean;
res.limitsAMH        = lims.audioMH;

% Dynamics of the Participant's Vocal Response
res.audioDynamics = niAn.audioDynamics;
end