function [sigs, trigs, vSigs] = dfMakePertSignal(trialLen, numTrial, sRateQ, sRateA, trialType, varargin)
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
%sRateQ:    Sampling Rate of device creating the stimulis (likely NIDAQ)
%sRateA:    Sampling Rate of program creating auditory cues (Audpater)
%trialType: Vector (length = numTrial) of order of trials (control/catch)
%varargin:  Any extra variables. Likely to just be a diagnostic flag

%sigs:  Per-trial digital signal outputted by NIDAQ
%trigs: Per-trial pertrubation start and stop points to be aligned with mic data

if isempty(varargin)
    diag = 0;
else
    diag = varargin{1};
end

minSt = 1.0; maxSt = 1.5;   %Hardset (1.0-1.5 seconds)
minLen = 1.0; maxLen = 1.5; %Hardset (1.0-1.5 seconds) 

trialLenP = trialLen*sRateQ; %Convert from seconds -> points for stimulus

valveBuffTime = 0.01; % seconds
valveBuffP    = valveBuffTime*sRateQ;

sigs     = zeros(trialLenP, numTrial);
trigs    = zeros(numTrial,2,3);
vSigs    = zeros(trialLenP, numTrial);
% Make a pert period for every trial, although only applied as a
% perturbatron signal, if actually a perturbed trial. See line 54
for i = 1:numTrial
    St_t   = (minSt + (maxSt-minSt)*rand);        % Start Time Seconds
    if diag == 1
        pLen_t = minLen;                          % Pert Len Seconds
    else
        pLen_t = (minLen + (maxLen-minLen)*rand); % Pert Len Seconds
    end
    Sp_t   = St_t + pLen_t;                       % Stop Time Seconds
    
    St_p   = round(sRateQ*St_t);                  % Start Time Points
    pLen_p = round(sRateQ*pLen_t);                % Pert Len Points  
    Sp_p   = St_p + pLen_p;                       % Stop TIme Points
    span   = St_p:Sp_p;                           % Start to Stop Span
    
    vSpan  = (St_p-valveBuffP):(Sp_p+valveBuffP); % Perturbation time with buffer for valve ON/OFF

    sig  = zeros(trialLenP,1);
    if trialType(i) ~= 0 % Only do this for perturbed trials
        sig(span) = 5;
    end
    
    vSig = zeros(trialLenP, 1);
    if trialType(i) == 2 % When TrialType is 2...open the valve, otherwise keep it in default position.
        vSig(vSpan) = 5;
    end
    
    sigs(:,i)    = sig;
    vSigs(:,i)   = vSig;
    trigs(i,:,1) = [St_t Sp_t];                  % Trigger in Seconds
    trigs(i,:,2) = [St_p Sp_p];                  % Trigger in Points (NIDAQ)
    trigs(i,:,3) = [St_p Sp_p]*(sRateA/sRateQ);  % Trigger in Points (Audapter)
end
end