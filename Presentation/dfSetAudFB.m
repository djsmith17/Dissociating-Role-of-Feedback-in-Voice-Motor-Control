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
    p.fb               = 1;
    p.bTimeDomainShift = 0;
    p.bCepsLift        = 0;
    p.dScale           = 1.00;
    p.nDelay           = 7;
    
    if isequal(lower(gender), 'female')
        p.pitchLowerBoundHz = 150;
        p.pitchUpperBoundHz = 300;
    elseif isequal(lower(gender), 'male')
        p.pitchLowerBoundHz = 80;
        p.pitchUpperBoundHz = 160;
    end 
    
    expParam.SSNw   = [];
    expParam.SSNfs  = [];
    
%PITCH-SHIFTED AUDITORY FEEDBACK OF THEIR VOICE
elseif expParam.AudFBSw == 1
    p.fb               = 1;
    p.bTimeDomainShift = 1;
    p.bCepsLift        = 1;
    p.dScale           = 1.00;
    p.nDelay           = 7; 
    
    if isequal(lower(gender), 'female')
        p.pitchLowerBoundHz = 150;
        p.pitchUpperBoundHz = 300;
    elseif isequal(lower(gender), 'male')
        p.pitchLowerBoundHz = 80;
        p.pitchUpperBoundHz = 160;
    end       
    
    expParam.SSNw   = [];
    expParam.SSNfs  = [];

%SPEECH-SHAPED MASKING NOISE
elseif expParam.AudFBSw == 2    
    p.fb               = 0;    % Audio File (Masking)
    p.bTimeDomainShift = 0;
    p.bCepsLift        = 0;
    p.dScale           = 1; % Headphone Scalar
    p.nDelay           = 7;
    
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