function [expParam, p] = dfSetAudFB(expParam, dirs, p)
%This function sets the type of Auditory Feedback to be played during an 
%experiment involving Audapter. If the Auditory feedback is speech-shaped
%masking noise, then you will need the masking noise wave file, SSN.wav.

%PITCH-SHIFTED AUDITORY FEEDBACK OF THEIR VOICE
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

    if fs ~= p.sr * p.downFact
        w = resample(w, p.sr * p.downFact, fs);              
    end
    if length(w) > maxPBSize
        w = w(1:maxPBSize);
    end
    Audapter('setParam', 'datapb', w, 1);
    
    expParam.SSNw   = w;
    expParam.SSNfs  = fs;
end 
end