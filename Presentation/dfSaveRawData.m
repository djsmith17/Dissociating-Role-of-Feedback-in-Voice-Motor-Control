function data = dfSaveRawData(expParam, dirs, p, audStimP, dataDAQ)
%Packages all of the experimental variables, and raw data into one 
%structure to be saved as a mat file and as wav files.

try
    data = AudapterIO('getData');
    
    data.expParam    = expParam; %Experimental Parameters
    data.dirs        = dirs;     %Directories
    data.p           = p;        %Audapter Parameters
    data.audStimP    = audStimP; %auditory stimulus Parameters
    data.DAQin       = dataDAQ;  %NIDAQ recordings ('Pressure/Force Sensors')
    save(fullfile(dirs.RecFileDir, [expParam.curSubCond dirs.saveFileSuffix]), 'data')

    audiowrite(fullfile(dirs.RecWaveDir,[expParam.curSubCond dirs.saveFileSuffix '_headOut.wav']), data.signalOut, expParam.sRateAnal)
    audiowrite(fullfile(dirs.RecWaveDir,[expParam.curSubCond dirs.saveFileSuffix '_micIn.wav']), data.signalIn, expParam.sRateAnal)
catch
    disp('Audapter decided not to show up today')
    data = [];
    return
end
end