function dfAnalysisJND()

JNDa.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
JNDa.participant  = 'Pilot10'; %List of multiple participants.

dirs = dfDirs(JNDa.project);
dirs.SavResultsDir = fullfile(dirs.Results, JNDa.participant, 'JND'); %Where to save results

if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

allRunData = [];
allmeanJND = [];
for ii = 1:4
    JNDa.run         = ['JNDpitch' num2str(ii)];
    
    dirs.SavFileDir  = fullfile(dirs.SavData, JNDa.participant, JNDa.run, ['ExperimentalParameters.mat']); %Where to find data
    
    load(dirs.SavFileDir)
    meanJND = dfAnalyzeThresholdJND(UD, 'reversals', 4)*.01; %Semitones
    
    allRunData = cat(1, allRunData, UD);
    allmeanJND = cat(1, allmeanJND, meanJND);
end

drawJNDResults(allRunData)

end