function data = dfSaveRawData(expParam, dirs)
%Packages all of the experimental variables, and raw data into one 
%structure to be saved as a mat file and as wav files.

try
    data = AudapterIO('getData');
    
    audiowrite(fullfile(dirs.RecWaveDir,[expParam.curExpTrial dirs.saveFileSuffix '_headOut.wav']), data.signalOut, expParam.sRateAnal)
    audiowrite(fullfile(dirs.RecWaveDir,[expParam.curExpTrial dirs.saveFileSuffix '_micIn.wav']), data.signalIn, expParam.sRateAnal)
catch
    disp('\nAudapter decided not to show up today')
    data = [];
    return
end
end