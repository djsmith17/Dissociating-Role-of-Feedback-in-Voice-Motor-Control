function dfAnalysisJND()
% dfAnalysisJND() loads a f0 JND data set (fAX#DRF.mat) and analyzed the raw
% data into interpretable results. After the results have been computed and
% organized, they are immediately plotted. Multiple JND results can be
% analyzed and plotted at the same time. By following the saved file
% structure of fAX1, fAX2,... etc, you can select the run numbers to
% analyze and view the run results next to each other in the plots.
%
% This function calls the following external functions
% -dfDirs
% -dfAnalyzeThresholdJND
% -drawJNDResults
%
% This script has the following subfunctions:
% -initJNDAnalysis()

close all

% Initalize the analysis structure
JNDa = initJNDAnalysis();

JNDa.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
JNDa.participants  = {'DRF1'}; % List of multiple participants.
JNDa.numPart       = length(JNDa.participants);
JNDa.runs          = {'fAX1'}; %List of multiple runs.
JNDa.numRuns       = length(JNDa.runs);

dirs = dfDirs(JNDa.project);
dirs.SavResultsDir = fullfile(dirs.Results, JNDa.participant, 'JND'); %Where to save results

if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

if strcmp(JNDa.runType, 'fAC')
    JNDa.tN = {'Diff'; 'Same'}; 
else
    JNDa.tN = {'First'; 'Last'};
end  

allJNDData  = [];
for ii = JNDa.runs2Analyze 
    curRun = [JNDa.runType num2str(ii)];
    JNDa.runs = cat(1, JNDa.runs, {curRun});
    
    dirs.SavFileDir  = fullfile(dirs.SavData, JNDa.participant, curRun, [JNDa.participant curRun 'DRF.mat']); %Where to find data
    
    fprintf('Loading Raw JND Data for %s %s\n', JNDa.participant, curRun)
    load(dirs.SavFileDir) % Returns UD
    
    JNDa.f0     = round(UD.subjf0, 1);
    JNDa.gender = UD.gender;
    
    JNDa.reversalsReached = cat(1, JNDa.reversalsReached, UD.reversals);
    JNDa.trialsCompleted  = cat(1, JNDa.trialsCompleted, UD.performedTrials);
    JNDa.timeElapsed      = cat(1, JNDa.timeElapsed, UD.elapsedTime);
    
    % Determine the JND Score and Accuracy of the last set of trials
    [JNDScore, lastSetAccu] = dfAnalyzeThresholdJND(UD, 'reversals', 4); %Cents
    
    JNDa.JNDScores       = cat(1, JNDa.JNDScores, JNDScore);
    JNDa.lastSetAccuracy = cat(1, JNDa.lastSetAccuracy, round(lastSetAccu, 1));
    JNDa.catchAccuracy   = cat(1, JNDa.catchAccuracy, round(UD.catchAccuracy));
    
    allJNDData  = cat(1, allJNDData, UD);
end

% Investigate some general stats now that we have compiled all the runs we
% want to look at.
JNDa = generalJNDStats(JNDa);

% Save the structure for future grouped analysis
dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [JNDa.participant JNDa.runType 'f0AcuityPooledResults.mat']);
fprintf('\nSaving Pooled JND Results for %s\n', JNDa.participant)
save(dirs.SavResultsFile, 'JNDa')

% Draw 
drawJNDResults(JNDa, dirs.SavResultsDir, allJNDData)
end

function JNDa = initJNDAnalysis()

JNDa.project     = [];
JNDa.participant = [];
JNDa.runType     = [];
JNDa.runs        = {};
JNDa.gender      = [];
JNDa.age         = [];
JNDa.f0          = [];

JNDa.instructions     = {};
JNDa.tN               = {};
JNDa.reversalsReached = [];
JNDa.trialsCompleted  = [];
JNDa.timeElapsed      = [];

JNDa.JNDScores       = [];
JNDa.lastSetAccuracy = [];
JNDa.catchAccuracy   = [];

JNDa.JNDScoreMean        = [];
JNDa.JNDScoreSE          = [];
JNDa.lastSetAccuracyMean = [];
JNDa.lastSetAccuracySE   = [];
end

function JNDa = generalJNDStats(JNDa)

JNDa.numJNDScores        = length(JNDa.JNDScores);

JNDa.JNDScoreMean        = round(mean(JNDa.JNDScores), 2);
JNDa.JNDScoreSE          = std(JNDa.JNDScores)/sqrt(JNDa.numJNDScores);
JNDa.lastSetAccuracyMean = round(mean(JNDa.lastSetAccuracy), 1);
JNDa.lastSetAccuracySE   = std(JNDa.lastSetAccuracy)/sqrt(JNDa.numJNDScores);
end