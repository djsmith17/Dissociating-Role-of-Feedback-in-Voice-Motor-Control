function rmsMean = dfCalcMeanRMS(rawData, varargin)
%Calculates the mean RMS based on Audapter recorded microphone samples

% Set a value for the baseline RMS. 
if isempty(varargin)
    rmsB  = 0.0000021689;
else
    rmsB = varargin{1};
end

rms     = rawData.rms(:,3);
rms     = rms(rms ~= 0);
rmsdB   = 20*log10(rms/rmsB);

%There were -Inf in my RMSdB. not exactly sure why, but this fixes the
%trial 1, lack of feeedback
rmsdB(rmsdB == -Inf) = 0;

rmsMean = mean(rmsdB);
end