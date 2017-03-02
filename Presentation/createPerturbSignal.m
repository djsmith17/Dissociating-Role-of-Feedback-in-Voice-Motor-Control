function [sigs, spans, spansT] = createPerturbSignal(trialLen, numTrial, sRateStim, trialType, expType)
%This function creates the digital signal for the NIDAQ needed to activate 
%the pertrubatron at each trial. It also keeps track of the time points 
%when the perturbation should activate and deactivate. This function 
%calculates these points for the sampling rate of the NIDAQ (s.Rate). 
%'spans' is then converted to the Audapter sampling rate (sRate) in the 
%main function.

%The perturbation starts at a semi-random time between 1.0s and 1.5s after 
%phonation. The pertrubation stops at a semi-random time between 
%1.0 and 1.5s after the start of pertrubation.

%trialLen:  The length of each trial in seconds
%numTrial:  The number of trials
%sRateStim: Sampling Rate of device creating the stimulis (likely NIDAQ)
%trialType: Vector (length = numTrial) of order of trials (control/catch)
%expType:   The string of the experiment that is being performed.

%sigs:  Per-trial digital signal outputted by NIDAQ
%spans: Per-trial pertrubation start and stop points to be aligned with mic data
%spansT: Per-trial pertrubation start and stop times to be aligned with mic data

expChk{1} = 'Somatosensory Perturbation_Perceptual';
expChk{2} = 'Auditory Perturbation_Perceptual';

minSt = 1.0; maxSt = 1.5;   %Hardset (1.0-1.5 seconds)
minLen = 1.0; maxLen = 1.5; %Hardset (1.0-1.5 seconds) 

trialLenP = trialLen*sRateStim; %Convert from seconds->points for stimulus

sigs    = zeros(trialLenP, numTrial);
spans   = zeros(numTrial,2);
spansT = zeros(numTrial,2);
for i = 1:numTrial
    St_t   = (minSt + (maxSt-minSt)*rand);    %Seconds
    pLen_t = (minLen + (maxLen-minLen)*rand); %Seconds
    Sp_t   = St_t + pLen_t;                   %Seconds
    
    St_p   = round(sRateStim*St_t);   %Points
    pLen_p = round(sRateStim*pLen_t); %Points  
    Sp_p   = St_p + pLen_p;        %Points
    span = St_p:Sp_p; 

    sig  = zeros(trialLenP,1);
    if trialType(i) == 1 && strcmp(expType, expChk{1}) %Only do this for SFPerturb
        sig(span) = 3; %For SFPerturb  (sometimes) =3. For AFPerturb always =0
    end
    
    sigs(:,i)    = sig;
    spans(i,:)   = [St_p Sp_p]; %Points
    spansT(i,:)  = [St_t Sp_t]; %Seconds
end
end