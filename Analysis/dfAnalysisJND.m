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
% -setCatchAcc
% -setDirection

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
    UD = setCatchAcc(UD);
    UD = setDirection(UD, JNDa, ii);
    
    JNDa.f0     = round(UD.subjf0, 1);
    JNDa.gender = UD.gender;
    
    JNDa.reversalsReached = cat(1, JNDa.reversalsReached, UD.reversals);
    JNDa.trialsCompleted  = cat(1, JNDa.trialsCompleted, UD.performedTrials);
    JNDa.timeElapsed      = cat(1, JNDa.timeElapsed, UD.elapsedTime);

    [meanJND, lastSetAccu] = dfAnalyzeThresholdJND(UD, 'reversals', 4); %Cents
    
    JNDa.JNDScores       = cat(1, JNDa.JNDScores, meanJND);
    JNDa.lastSetAccuracy = cat(1, JNDa.lastSetAccuracy, round(lastSetAccu, 1));
    JNDa.catchAccuracy   = cat(1, JNDa.catchAccuracy, round(UD.catchAccuracy));
    
    allJNDData  = cat(1, allJNDData, UD);
end

JNDa.JNDScoreMean        = round(mean(JNDa.JNDScores), 2);
JNDa.lastSetAccuracyMean = round(mean(JNDa.lastSetAccuracy), 1);

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

function UD = setCatchAcc(UD)
% setCatchAcc() calculates the stats on performed catch trials, only if the
% JND task had catch trials. Not all versions do.

if ~isfield(UD, 'catchTrials')
    UD.performedTrials = length(UD.catchResponse);
    UD.JNDTrials       = length(UD.reversal);
    UD.catchTrials     = length(UD.catchResponse) - length(UD.reversal);
    UD.reversals       = max(UD.reversal);
    UD.catchCorrect    = sum(UD.catchResponse == 0);
    UD.catchAccuracy   = 100*(UD.catchCorrect/UD.catchTrials);
end
end

function UD = setDirection(UD, JNDa, ii)
% setDirection is a function I made because I goofed in recording some
% pilot data. You can ignore this otherwise. 

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