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

prompt = {'Subject ID:',...
          'Run Type:',...
          'Runs to Analyze:'};
name = 'Subject Information';
numlines = 1;
defaultanswer = {'null', 'fAX', '1, 2, 3, 4'};
answer = inputdlg(prompt, name, numlines, defaultanswer);

if isempty(answer)
    return
end

% Initalize the analysis structure
JNDa = initJNDAnalysis();

JNDa.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
JNDa.participant  = answer{1};
JNDa.runType      = answer{2};
num = textscan(answer{3}, '%f', 'Delimiter', ',');
JNDa.runs2Analyze = num{1}';

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

JNDa.JNDScoreMean        = round(mean(JNDa.JNDScores), 2);
JNDa.lastSetAccuracyMean = round(mean(JNDa.lastSetAccuracy), 1);

dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [JNDa.participant JNDa.runType 'f0AcuityPooledResults.mat']);
fprintf('\nSaving Pooled JND Results for %s\n', JNDa.participant)
save(dirs.SavResultsFile, 'JNDa')

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
JNDa.reversalsReached = [];
JNDa.trialsCompleted  = [];
JNDa.timeElapsed      = [];

JNDa.JNDScores       = [];
JNDa.lastSetAccuracy = [];
JNDa.catchAccuracy   = [];

JNDa.JNDScoreMean        = [];
JNDa.JNDScoreSD          = [];
JNDa.lastSetAccuracyMean = [];
JNDa.lastSetAccuracySD   = [];
end