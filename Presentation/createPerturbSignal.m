function [sigs, spans] = createPerturbSignal(s, numTrial, trialLen, trialType)
%This function creates the digital signal for the NIDAQ needed to activate 
%the pertrubatron at each trial. It also keeps track of the time points 
%when the perturbation should activate and deactivate. This function 
%calculates these points for the sampling rate of the NIDAQ (s.Rate). 
%'spans' is then converted to the Audapter sampling rate (sRate) in the 
%main function.

%The perturbation starts at a semi-random time between 1.7s and 2.1s after 
%phonation. The pertrubation stops at a semi-random time between 
%1.0 and 1.5s after the start of pertrubation.

%s:         NIDAQ object handle for keeping track of variables and I/O
%numTrial:  The number of trials
%trialLen:  The length of each trial in points
%trialType: Vector (length = numTrial) of order of trials (control/catch)

%sigs:  Per-trial digital signal outputted by NIDAQ
%spans: Per-trial pertrubation start and stop points to be aligned with mic data

minSt = 1.7; maxSt = 2.1;   %Hardset (1.7-2.1 seconds)
minLen = 1.0; maxLen = 1.5; %Hardset (1.0-1.5 seconds) 

sigs  = zeros(trialLen, numTrial);
spans = zeros(numTrial,2);
for i = 1:numTrial
    St   = round(s.Rate*(minSt + (maxSt-minSt)*rand)); 
    pLen = round(s.Rate*(minLen + (maxLen-minLen)*rand));   
    Sp   = St+pLen;    
    span = St:Sp; 

    sig  = zeros(trialLen,1);
    if trialType(i) == 1
        sig(span) = 3; %For SFPerturb =3. For AFPerturb =0
    end
    
    sigs(:,i)  = sig;
    spans(i,:) = [St Sp];
end
end