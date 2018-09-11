function [p, SSNw, SSNfs] = dfSetAudFB(expParam, dirs, p)
% [expParam, p] = dfSetAudFB(expParam, dirs, p) sets the type of Auditory 
% feedback to be used in experiments investigating voice motor control 
% using Audapter. 
% 
% This updated the expParam and p structures and then passes them back out
% to be used in the main experiment. 
%
% This handles three different cases of Audiotory feedback
% NORMAL AUDITORY FEEDBACK OF THEIR VOICE
% PITCH-SHIFTED AUDITORY FEEDBACK OF THEIR VOICE
% SPEECH-SHAPED MASKING NOISE
%
% This function has a subfunction below named: setLoudRatio
  
dB           = expParam.headGain;
gender       = expParam.gender;
p.bBypassFmt = 1;               % No Formant tracking
p.dScale     = setLoudRatio(dB);% Scale of output from input
% p.nDelay     = 7;

if isequal(lower(gender), 'female')
    p.pitchLowerBoundHz = 150;
    p.pitchUpperBoundHz = 300;
elseif isequal(lower(gender), 'male')
    p.pitchLowerBoundHz = 80;
    p.pitchUpperBoundHz = 160;
end 

SSNw   = [];
SSNfs  = [];

%NORMAL AUDITORY FEEDBACK OF THEIR VOICE
if  expParam.AudFBSw == 0
    p.fb                = 1; % Microphone FB
    p.bTimeDomainShift  = 0; % No pitch-shifting
    p.bCepsLift         = 0;
    
%PITCH-SHIFTED AUDITORY FEEDBACK OF THEIR VOICE
elseif expParam.AudFBSw == 1
    p.fb                = 1; % Microphone FB
    p.bTimeDomainShift  = 1; % Pitch-shifting
    p.bCepsLift         = 1;

%SPEECH-SHAPED MASKING NOISE
elseif expParam.AudFBSw == 2
    p.fb                = 0; % Audio File (Masking)
    p.bTimeDomainShift  = 0; % No pitch-shifting
    p.bCepsLift         = 0;
    
    % For this trial (or set of trials), how long do we need noise?
    noiseTime = calcMaskLen(expParam);
    
    % Generate a full length masking noise signal for the length we need
    [w, fs] = createSessionNoise(dirs, noiseTime);
    
%     [w, fs] = audapterGeneratedNoise(dirs, p);

    SSNw   = w;
    SSNfs  = fs;
else
    error('ERROR in dfSetAudFB: Inappropriate feedback method selected')
end 
end

function dScale = setLoudRatio(dB)
% dB should be a positive or negative decimal value of change in dB
% representing how much you want to scale the headphones against the input 
% microphone level.

dScale = 10^(dB/20);
end

function [w, fs] = audapterGeneratedNoise(dirs, p)

%Uses Speech-Shaped Noise stored in util
noiseWavFN = fullfile(dirs.Prelim, 'SSN.wav'); 

maxPBSize  = Audapter('getMaxPBLen');

check_file(noiseWavFN);
[w, fs] = read_audio(noiseWavFN);

if fs ~= p.sr * p.downFact
    w = resample(w, p.sr * p.downFact, fs);              
end
if length(w) > maxPBSize
    w = w(1:maxPBSize);
end
Audapter('setParam', 'datapb', w, 1);
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
noiseLen = round(noiseTime*fs);

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