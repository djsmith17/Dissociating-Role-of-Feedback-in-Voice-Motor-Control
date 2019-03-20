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
% -AnalyzeRawJNDData()
tic
close all

% Initalize the analysis structure
JNDa.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
JNDa.participants = {'DRF1',...
                      'DRF2',...
                      'DRF4',...
                      'DRF5',...
                      'DRF6',...
                      'DRF7',...
                      'DRF8',...
                      'DRF9',...
                      'DRF10',...
                      'DRF12',...
                      'DRF13',...
                      'DRF14',...
                      'DRF15',...
                      'DRF16',...
                      'DRF17',...
                      'DRF18',...
                      'DRF19'}; % List of multiple participants.
JNDa.numPart      = length(JNDa.participants);
JNDa.runs         = {'fAX1', 'fAX2','fAX3','fAX4'}; %List of multiple runs.
JNDa.numRuns      = length(JNDa.runs);

dirs = dfDirs(JNDa.project);

for jj = 1:JNDa.numPart
    curPart = JNDa.participants{jj}; % Current Participant
    dirs.SavResultsDir = fullfile(dirs.Results, curPart, 'JND'); %Where to save results

    if exist(dirs.SavResultsDir, 'dir') == 0
        mkdir(dirs.SavResultsDir) %If the folder we are saving does not exist, let's make it
    end

    for ii = 1:JNDa.numRuns
        curRun     = JNDa.runs{ii}; % Current Run
        dirs.SavFileDir = fullfile(dirs.SavData, curPart, curRun, [curPart curRun 'DRF.mat']); %Where to find raw data

        fprintf('****%s  %s****\n', curPart, curRun)
        fprintf('Loading Raw JND Data\n')
        load(dirs.SavFileDir) % Returns UD
        
        % Organize the raw data into a result structure named 'resJND'
        resJND = AnalyzeRawJNDData(UD);
        dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [curPart curRun 'ResultsDRF.mat']); %What to name the organized results
        fprintf('Saving Individual JND Results\n\n')
        save(dirs.SavResultsFile, 'resJND')    
    end
end
fprintf('Elapsed time was %f min\n', toc/60)
end

function resJND = initJNDAnalysis()

resJND.participant = [];
resJND.gender      = [];
resJND.f0          = [];
resJND.run         = [];

resJND.instructions     = {};
resJND.selectOpt        = {};

resJND.reversalsReached = [];
resJND.trialsCompleted  = [];
resJND.timeElapsed      = [];

resJND.distProgression = [];

resJND.trialsAtReversals = [];
resJND.distAtReversals   = [];

resJND.trialsAtCorrectOpt1   = [];
resJND.distAtCorrectOpt1     = [];
resJND.trialsAtIncorrectOpt1 = [];
resJND.distAtIncorrectOpt1   = [];
resJND.trialsAtCorrectOpt2   = [];
resJND.distAtCorrectOpt2     = [];
resJND.trialsAtIncorrectOpt2 = [];
resJND.distAtIncorrectOpt2   = [];

resJND.JNDScore        = [];
resJND.LastSetAccuracy = [];
resJND.catchAccuracy   = [];
end

function resJND = AnalyzeRawJNDData(UD)
% This calls the following subfunctions
% -initJNDAnalysis()
% -dfAnalyseThresholdJND()

fprintf('Starting Individual JND Analysis\n')

resJND = initJNDAnalysis();

resJND.participant = UD.subject;
resJND.gender      = UD.gender;
resJND.f0          = round(UD.subjf0, 1);
resJND.run         = UD.run;
    
resJND.instructions = UD.inst;    
resJND.selectOpt = {'First'; 'Last'};

resJND.reversalsReached = UD.reversals;
resJND.trialsCompleted  = UD.performedTrials;
resJND.timeElapsed      = UD.elapsedTime;

resJND.distProgression = UD.x;

resJND.trialsAtReversals = find(UD.reversal~=0);
resJND.distAtReversals   = UD.x(resJND.trialsAtReversals);

resJND.trialsAtCorrectOpt1   = find(UD.allTrialTypes == 1);
resJND.trialsAtIncorrectOpt1 = find(UD.allTrialTypes == 2);
resJND.trialsAtCorrectOpt2   = find(UD.allTrialTypes == 3);
resJND.trialsAtIncorrectOpt2 = find(UD.allTrialTypes == 4);

resJND.distAtCorrectOpt1   = UD.x(resJND.trialsAtCorrectOpt1)';
resJND.distAtIncorrectOpt1 = UD.x(resJND.trialsAtIncorrectOpt1)';
resJND.distAtCorrectOpt2   = UD.x(resJND.trialsAtCorrectOpt2)';
resJND.distAtIncorrectOpt2 = UD.x(resJND.trialsAtIncorrectOpt2)';

% Determine the JND Score and Accuracy of the last set of trials
[JNDScore, LastSetAccuracy] = dfAnalyzeThresholdJND(UD, 'reversals', 4); %Cents
resJND.JNDScore        = JNDScore;
resJND.LastSetAccuracy = round(LastSetAccuracy, 1);
resJND.catchAccuracy   = round(0); %Currently N/A, but maybe used again later
end