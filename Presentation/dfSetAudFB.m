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

gender = expParam.gender;

%NORMAL AUDITORY FEEDBACK OF THEIR VOICE
if  expParam.AudFBSw == 0
    p.fb          = 1;  % Microphone FB
    p.bPitchShift = 0;  % No pitch-shifting
    p.bBypassFmt  = 1;  % No Formant tracking
    p.dScale      = 1.00;
    
    expParam.SSNw   = [];
    expParam.SSNfs  = [];
    
%PITCH-SHIFTED AUDITORY FEEDBACK OF THEIR VOICE
elseif expParam.AudFBSw == 1
    p.bTimeDomainShift = 1;
    
    if isequal(lower(gender), 'female')
        p.pitchLowerBoundHz = 150;
        p.pitchUpperBoundHz = 300;
    elseif isequal(lower(gender), 'male')
        p.pitchLowerBoundHz = 80;
        p.pitchUpperBoundHz = 160;
    end
    
    p.frameLen = 64;
    p.nDelay = 7;
    p.bCepsLift = 1;
    
%     p.timeDomainPitchShiftSchedule = [0, 1.0; 1, 1.0; 1.150, 1.0; 2.0, 1.0; 2.15, 1.0];
    p.timeDomainPitchShiftSchedule = [0, 1.0; 1, 1.0; 1.150, 0.9438; 2.0, 0.9438; 2.15, 1.0];
    p.rmsThresh = 0.011;
    
%     p.fb          = 1;    % Microphone FB
%     p.bPitchShift = 1;    % Pitch-shifting
%     p.bBypassFmt  = 1;    % No Formant tracking
%     p.dScale      = 1.00; %Headphone Scalar
    
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