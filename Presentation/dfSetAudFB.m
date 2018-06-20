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
    p.bBypassFmt  = 1;  % No Formant tracking
    p.dScale      = 2.24;
    
    expParam.SSNw   = [];
    expParam.SSNfs  = [];
    
%PITCH-SHIFTED AUDITORY FEEDBACK OF THEIR VOICE
elseif expParam.AudFBSw == 1
    p.fb          = 1;    % Microphone FB
    p.bPitchShift = 1;    % Pitch-shifting
    p.bBypassFmt  = 1;    % No Formant tracking
    p.dScale      = 2.24; %Headphone Scalar
    
    expParam.SSNw   = [];
    expParam.SSNfs  = [];

%SPEECH-SHAPED MASKING NOISE
elseif expParam.AudFBSw == 2
    p.fb          = 0;    % Audio File (Masking)
    p.bPitchShift = 0;
    p.bBypassFmt  = 1;    % No Formant tracking
    p.dScale      = 1; % Headphone Scalar
    
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
    
    expParam.SSNw   = w;
    expParam.SSNfs  = fs;
else
    disp('ERROR in dfSetAudFB: Inappropriate feedback method selected')
    return
end 
end