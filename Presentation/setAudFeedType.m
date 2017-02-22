function p = setAudFeedType(p, masking)
%This function sets the type of Auditory Feedback to be played during an 
%experiment involving Audapter. If the Auditory feedback is speech-shaped
%masking noise, then you will need the masking noise wave file, SSN.wav.

p.masking = masking;

if masking == 0
    p.fb          = 1;
    p.bPitchShift = 1;
    p.dScale      = 1; %Headphone Scalar
elseif masking == 1
    p.fb          = 2;
%     p.fb3Gain     = 2.0;
    p.bPitchShift = 0;
    p.dScale      = 1; %Headphone Scalar
    noiseWavFN = 'util\SSN.wav'; %Uses Speech-Shaped Noise stored in util
    
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
end