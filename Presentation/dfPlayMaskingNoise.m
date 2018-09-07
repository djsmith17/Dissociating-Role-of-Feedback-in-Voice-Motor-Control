function dfPlayMaskingNoise()
%Test the levels of masking noise based on how the paradigm is set up for
%voice perturbations
%
% This script calls the following 8 functions:
% dfDirs.m
% dfSetAudFB.m
%
% This uses the toolbox from MATLAB-Toolboxes
% speechres

recordedRMS = (5+rand(1991,1))/1000;
rmsMean     = mean(recordedRMS);

expParam = dfInitExpParam();

expParam.curSess = 'Masking Noise Test';

expParam.numTrial   = 2;
expParam.AudFB      = 'AC Masking Noise';
expParam.AudFBSw    = 2; % Masking Noise

% Set our dirs based on the project
dirs = dfDirs(expParam.project);

% Set up Audapter Parameters and function paths
[expParam, p] = dfInitAudapater(dirs, expParam);

% Set up Auditory Feedback (Masking Noise)
[expParam, p] = dfSetAudFB(expParam, dirs, p);

% This is where the fun begins
fprintf('\nStarting Trials\n\n')

% Dim the lights (Set the visual Feedback)
[msrStr, annoStr] = dfSetVisFB(expParam.curSess, expParam.targRMS, expParam.boundsRMS);

noiseTime = calcMaskLen(expParam);
[sessionNoise, noiseFs] = createSessionNoise(dirs, noiseTime);

sound(sessionNoise, noiseFs)

ET = tic;
% Open the curtains
pause(expParam.rdyPause); % Let them breathe a sec
set(annoStr.Ready, 'Visible','off'); % Turn off 'Ready?'
for ii = 1:expParam.numTrial

    
    %Set the OST and PCF functions
    Audapter('ost', expParam.ostFN, 0);
    Audapter('pcf', expParam.pcfFN, 0);
    
    %Cue to begin trial
    set(annoStr.plus, 'Visible','on');
    pause(expParam.cuePause)
    set(annoStr.plus, 'Visible','off');
    
    %Phonation Start
    set([annoStr.EEE annoStr.visTrig],'Visible','on');
    
    fprintf('Trial %d\n', ii)
    AudapterIO('init', p);
    Audapter('reset');
    Audapter('start');
    pause(expParam.buffPause)

    pause(expParam.trialLen)
    
    %Phonation End
    set([annoStr.EEE annoStr.visTrig],'Visible','off');
    pause(expParam.endPause)
    Audapter('stop');
    
    % Compare against baseline and updated Visual Feedback
    [color, newPos, loudResult] = dfUpdateVisFB(msrStr, rmsMean);
     
    % Provide Loudness Feedback
    set(annoStr.LoudRec, 'position', newPos);
    set(annoStr.LoudRec, 'Color', color, 'FaceColor', color);
    set([annoStr.LoudRec annoStr.fbLines], 'Visible', 'on');    
    pause(expParam.resPause)
    set([annoStr.LoudRec annoStr.fbLines], 'Visible', 'off');
end
close all

elapsed_time = toc(ET);   % Elapsed Time of the experiment
fprintf('\nElapsed Time: %f (s)\n', elapsed_time)
end

function noiseTime = calcMaskLen(expParam)

numTrial = expParam.numTrial;

rdyTime  = expParam.rdyPause;  % Ready Message
cueTime  = expParam.cuePause;  % Cue period
buffTime = expParam.buffPause; % Buffer to begin phonating
trlTime  = expParam.trialLen;  % Phonation period
endTime  = expParam.endPause;  % Buffer to end phonating
resTime  = expParam.resPause;  % Rest/Feedback period

noiseTime = rdyTime + (cueTime + buffTime + trlTime + endTime + resTime)*numTrial + 2;
end

function [sessionNoise, fs] = createSessionNoise(dirs, noiseTime)

maskFile = fullfile(dirs.Prelim, 'SSN.wav');

[wavFile, fs] = audioread(maskFile);
wavLen   = length(wavFile);
noiseLen = noiseTime*fs;

rampUpSp = round(2*fs) + 1;
rampDnSt = round((noiseTime-2)*fs);

rampUpIdx = 1:rampUpSp;
rampUpL   = length(rampUpIdx);
rampDnIdx = rampDnSt:noiseLen;
rampDnL   = length(rampDnIdx);

rampUp = linspace(0, 1, rampUpL);
rampDn = linspace(1, 0, rampDnL);

numRep   = noiseLen/wavLen; % How many repetitions of the .wav file (decimal)
minInt   = floor(numRep);   % Min number of whole repeitions (integer)
noiseInt = repmat(wavFile', [1, minInt]);

remRep   = numRep - minInt; % How many decimal amounts left?
remIdx   = round(wavLen*remRep);
noiseRem = wavFile(1:remIdx)';

fullNoise = [noiseInt noiseRem];

rampFilt = ones(size(fullNoise));
rampFilt(rampUpIdx) = rampUp;
rampFilt(rampDnIdx) = rampDn;

sessionNoise = fullNoise.*rampFilt;
end