function dBSPLTrialMean = dfCalcMeanRMS(rawData, varargin)
% Calculates the mean dB SPL from the rms of a recorded trial

% Reference Sound Pressure
if isempty(varargin)
    refSPress  = 0.00002;    % Default 20 micropascals for air
else
    refSPress = varargin{1};
end

rmsVec     = rawData.rms(:,1);           % smoothed rms from Audapter
dBSPLVec   = 20*log10(rmsVec/refSPress); % convert to dB

% Remove any -Inf in the vector
dBSPLVec(dBSPLVec == -Inf) = 0;

dBSPLTrialMean = mean(dBSPLVec); % Average of the vector
end