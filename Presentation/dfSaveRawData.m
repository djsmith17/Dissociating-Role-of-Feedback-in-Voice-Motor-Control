function data = dfSaveRawData(expParam, dirs)
% data = dfSaveRawData(expParam, dirs) loads all the Audapter data that is 
% currently sitting in the buffer and makes it available to be processed
% This function also immediately saves .wav files of the microphone and 
% headphone data from the raw Audapter data.
%
% This function requires the current dirs (dirs) for saving Wav Files and 
% the experimental parameters (expParam) associated with the current trial.

fs        = expParam.sRateAnal;
sessTrial = expParam.curSessTrial;

try
    data = AudapterIO('getData');
    
    audiowrite(fullfile(dirs.RecWaveDir,[sessTrial dirs.saveFileSuffix '_headOut.wav']), data.signalOut, fs)
    audiowrite(fullfile(dirs.RecWaveDir,[sessTrial dirs.saveFileSuffix '_micIn.wav']), data.signalIn, fs)
catch
    fprintf('\nAudapter decided not to show up today')
    data = [];
    return
end
end