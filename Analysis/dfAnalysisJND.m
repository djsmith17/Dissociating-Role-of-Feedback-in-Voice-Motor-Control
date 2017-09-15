function dfAnalysisJND()

JNDa.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
JNDa.participant  = 'Pilot11'; %List of multiple participants.
runs2Analyze      = 1:2;

dirs = dfDirs(JNDa.project);
dirs.SavResultsDir = fullfile(dirs.Results, JNDa.participant, 'JND'); %Where to save results

if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

allRunData = [];
allmeanJND = [];
allCatchAcc = [];
for ii = runs2Analyze 
    JNDa.run         = ['fAC' num2str(ii)];
    
    dirs.SavFileDir  = fullfile(dirs.SavData, JNDa.participant, JNDa.run, [JNDa.participant JNDa.run 'DRF.mat']); %Where to find data
    
    load(dirs.SavFileDir)
    UD = setCatchAcc(UD);
    UD = setDirection(UD, JNDa, ii);
    meanJND = dfAnalyzeThresholdJND(UD, 'reversals', 4); %Cents
    
    allRunData = cat(1, allRunData, UD);
    allmeanJND = cat(1, allmeanJND, meanJND);
    allCatchAcc = cat(1, allCatchAcc, UD.catchAccuracy);
end

drawJNDResults(JNDa, dirs, runs2Analyze, allRunData, allmeanJND, allCatchAcc)
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

function UD = setDirection(UD, JNDa, ii)
if ~isfield(UD, 'JNDDirection')
    if strcmp(JNDa.participant, 'Pilot9')
        if ii == 1 || ii == 3
            UD.JNDDirection = 'Above';
        else
            UD.JNDDirection = 'Below';
        end
    elseif strcmp(JNDa.participant, 'Pilot10')
        if ii == 2 || ii == 4
            UD.JNDDirection = 'Above';
        else
            UD.JNDDirection = 'Below';
        end 
    end
end
end