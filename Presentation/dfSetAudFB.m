function [expParam, p] = dfSetAudFB(expParam, dirs, p)
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

%NORMAL AUDITORY FEEDBACK OF THEIR VOICE
if  expParam.AudFBSw == 0
    p.fb          = 1;  % Microphone FB
    p.bPitchShift = 0;  % No pitch-shifting
    
    expParam.SSNw   = [];
    expParam.SSNfs  = [];
    
%PITCH-SHIFTED AUDITORY FEEDBACK OF THEIR VOICE
elseif expParam.AudFBSw == 1
    p.fb          = 1;   % Microphone FB
    p.bPitchShift = 1;   % Pitch-shifting
%     p.dScale      = 1; %Headphone Scalar
    
    expParam.SSNw   = [];
    expParam.SSNfs  = [];

%SPEECH-SHAPED MASKING NOISE
elseif expParam.AudFBSw == 2
    p.fb          = 2;   % Audio File (Masking)
    p.bPitchShift = 0;
%     p.dScale      = 1; %Headphone Scalar
    
    %Uses Speech-Shaped Noise stored in util
    noiseWavFN = fullfile(dirs.Prelim, 'SSN.wav'); 
    
    maxPBSize  = Audapter('getMaxPBLen');

    check_file(noiseWavFN);
    [w, fs] = read_audio(noiseWavFN);
    
    rampLen  = expParam.AFRampLen;
    rampLenP = rampLen*fs;
    ramp = linspace(0, 1, rampLenP);
    win = ones(size(w));
    win(1:rampLenP) = ramp;
    wWin = w.*win;

    if fs ~= p.sr * p.downFact
        wWin = resample(wWin, p.sr * p.downFact, fs);              
    end
    if length(wWin) > maxPBSize
        wWin = wWin(1:maxPBSize);
    end
    Audapter('setParam', 'datapb', wWin, 1);
    
    expParam.SSNw   = wWin;
    expParam.SSNfs  = fs;
else
    disp('ERROR in dfSetAudFB: Inappropriate feedback method selected')
    return
end 
end