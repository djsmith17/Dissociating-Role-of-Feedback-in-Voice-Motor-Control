function rmsMean = dfCalcMeanRMS(rawData, varargin)
%Calculates the mean RMS based on Audapter recorded microphone samples

if isempty(varargin)
    refSPL  = 0.00002;  % 20 micropascals for air
else
    refSPL = varargin{1};
end

rms     = rawData.rms(:,3);
rmsdB   = 20*log10(rms/refSPL);

%There were -Inf in my RMSdB. not exactly sure why, but this fixes the
%trial 1, lack of feeedback
rmsdB(rmsdB == -Inf) = 0;

rmsMean = mean(rmsdB);
end