function dfAnalysisJND()

JNDa.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
JNDa.participant  = 'Pilot9'; %List of multiple participants.

dirs = dfDirs(JNDa.project);
dirs.SavResultsDir = fullfile(dirs.Results, JNDa.participant, 'JND'); %Where to save results

if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

allRunData = [];
allmeanJND = [];
allCatchAcc = [];
for ii = 1:4
    JNDa.run         = ['fA' num2str(ii)];
    
    dirs.SavFileDir  = fullfile(dirs.SavData, JNDa.participant, JNDa.run, [JNDa.participant JNDa.run 'DRF.mat']); %Where to find data
    
    load(dirs.SavFileDir)
    UD = setCatchAcc(UD);
    meanJND = dfAnalyzeThresholdJND(UD, 'reversals', 4); %Cents
    
    allRunData = cat(1, allRunData, UD);
    allmeanJND = cat(1, allmeanJND, meanJND);
    allCatchAcc = cat(1, allCatchAcc, UD.catchAccuracy);
end

drawJNDResults(JNDa, dirs, allRunData, allmeanJND)
end

function UD = setCatchAcc(UD)
if ~isfield(UD, 'catchTrials')
    UD.performedTrials = length(UD.catchResponse);
    UD.JNDTrials = length(UD.reversal);
    UD.catchTrials = length(UD.catchResponse) - length(UD.reversal);
    UD.reversals = max(UD.reversal);
    UD.catchCorrect = sum(UD.catchResponse == 0);
    UD.catchAccuracy = 100*(UD.catchCorrect/UD.catchTrials);
end
end