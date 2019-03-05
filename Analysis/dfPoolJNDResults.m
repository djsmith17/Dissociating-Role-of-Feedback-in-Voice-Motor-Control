function dfPoolJNDResults()

close all
pA.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control'; 
pA.pAnalysis     = 'DRF_JND'; % Change this name to load different pooled data sets Ex: SfN2017, LarynxPos

dirs               = dfDirs(pA.project);
dirs.SavResultsDir = fullfile(dirs.Results, 'Pooled Analyses', pA.pAnalysis);        % The directory to save Files/Figures
dirs.PooledConfigF = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'PooledConfig.mat']);% The Pooled Config File for this Pooled Analysis

% Can we find the Pooled Config File?
if exist(dirs.PooledConfigF, 'file') == 0
    fprintf('\nERROR: Pooled Config File %s does not exist! Please create it with a GenConfig Function\n', dirs.PooledConfigF)
    return
else
    % Load the configuration file. Should return a data structure 'cF'
    load(dirs.PooledConfigF)
end 

pA.participants  = cF.participants; % List of multiple participants.
pA.numPart       = length(pA.participants);
pA.runs          = cF.runs;         % All runs to consider 
pA.numRuns       = length(pA.runs);
pA.cond          = cF.cond;         % Conditions to test against
pA.numCond       = length(pA.cond); 
pA.condVar       = cF.condVar;      % Variable to test the condition
pA.testExt       = cF.testExt;

pA.pltNameMVi    = cell(pA.numPart, 1);
pA.pltNameMVm    = [pA.pAnalysis 'MeanSubj' pA.testExt];

% Load all saved results and order into a large data structure
allDataStr = []; % (numPart x numRuns)
for ii = 1:pA.numPart
    participant = pA.participants{ii};
    fprintf('Loading Runs for %s\n', participant)
    
    pA.pltNameMVi{ii} = [pA.pAnalysis participant pA.testExt];
    subjRes  = [];
    for jj = 1:pA.numRuns
        run              = pA.runs{jj};
        dirs.SavFileDir  = fullfile(dirs.Results, participant, 'JND');                      % Where results are saved
        dirs.SavFile     = fullfile(dirs.SavFileDir, [participant run 'ResultsDRF.mat']); % Run Results file to load

        if exist(dirs.SavFile, 'file') == 0
            error('ERROR: The Results file for Run %s %s does not exist yet\n', participant, run)
        else   
            load(dirs.SavFile)
            % Returns a results struture of 'res'
        end
        
        subjRes = cat(2, subjRes, resJND);
    end
    allDataStr = cat(1, allDataStr, subjRes);
end

allSubjRes = initSortedStruct(pA.numCond, pA.numRuns);
allSubjRes.subject = 'Mean Participant Response';
allSubjRes.curSess = allSubjRes.subject;
allSubjRes.cond    = pA.cond;

for ii = 1:pA.numPart
    participant = pA.participants{ii};
    fprintf('Sorting task conditions for %s\n', participant)
    
    sortStruc         = initSortedStruct(pA.numCond, pA.numRuns);
    sortStruc.subject = ['Participant ' num2str(ii)]; % Pooled Analysis Name
    
    sortStruc.curSess = sortStruc.subject;
 
    for jj = 1:pA.numRuns
        curRes = allDataStr(ii, jj);

        sortStruc.studyID = curRes.participant; % Study ID
        sortStruc.gender  = curRes.gender;
%         sortStruc.age     = curRes.age;
        sortStruc.f0      = curRes.f0;
        
        sortStruc  = combineCondTrials(pA, curRes, sortStruc);       
        allSubjRes = combineCondTrials(pA, curRes, allSubjRes);
        
        sortStruc.distProgression{jj}   = curRes.distProgression;
        sortStruc.trialsAtReversals{jj} = curRes.trialsAtReversals;
        sortStruc.distAtReversals{jj}   = curRes.distAtReversals;
        
        sortStruc.trialsAtCorrectOpt1{jj}   = curRes.trialsAtCorrectOpt1;
        sortStruc.distAtCorrectOpt1{jj}     = curRes.distAtCorrectOpt1;
        sortStruc.trialsAtIncorrectOpt1{jj} = curRes.trialsAtIncorrectOpt1;
        sortStruc.distAtIncorrectOpt1{jj}   = curRes.distAtIncorrectOpt1;
        sortStruc.trialsAtCorrectOpt2{jj}   = curRes.trialsAtCorrectOpt2;
        sortStruc.distAtCorrectOpt2{jj}     = curRes.distAtCorrectOpt2;
        sortStruc.trialsAtIncorrectOpt2{jj} = curRes.trialsAtIncorrectOpt2;
        sortStruc.distAtIncorrectOpt2{jj}   = curRes.distAtIncorrectOpt2;
    end
        
    sortStruc = meanCondTrials(sortStruc);
    sortStruc.pltName = pA.pltNameMVi{ii};
    
    pooledRunStr(ii)   = sortStruc;
    
    allSubjRes.expType    = sortStruc.expType;
    allSubjRes.gender{ii} = sortStruc.gender;
%     allSubjRes.age(ii)    = sortStruc.age;
    
%     Save the structure for future grouped analysis
    dirs.SavResultsDirParti = fullfile(dirs.Results, participant, 'JND');
    dirs.SavResultsFile = fullfile(dirs.SavResultsDirParti, [participant 'f0AcuityPooledResults.mat']);
    fprintf('\nSaving Pooled JND Results for %s\n', participant)
    save(dirs.SavResultsFile, 'sortStruc')

    % Draw 
%     drawJNDResults(sortStruc, dirs.SavResultsDirParti)
    close all
end

stat.participants = {pooledRunStr.subject}';
stat.genders      = {pooledRunStr.gender}';
stat.f0           = {pooledRunStr.f0}';
stat.JNDScoreMean = {pooledRunStr.JNDScoreMean}';
stat.JNDScoreSE   = {pooledRunStr.JNDScoreSE}';
stat.lastSetAccuracyMean = {pooledRunStr.lastSetAccuracyMean}';
stat.lastSetAccuracySE   = {pooledRunStr.lastSetAccuracySE}';

JNDStatTable = table(stat.participants,...
                     stat.genders,...
                     stat.f0,...
                     stat.JNDScoreMean,...
                     stat.JNDScoreSE,...
                     stat.lastSetAccuracyMean,...
                     stat.lastSetAccuracySE,...
                     'VariableNames',...
                     {'Participant', 'gender', 'f0', 'JNDScoreMean', 'JNDScoreSE', 'lastSetAccuracyMean', 'lastSetAccuracySE'});
                 
figure
hist(cell2mat(JNDStatTable.JNDScoreMean))
title('Mean JND Score Distribution')

figure
hist(cell2mat(JNDStatTable.lastSetAccuracyMean))
title('Mean Last 4 Reversals Accuracy Distribution')

JNDTableCSV = fullfile(dirs.SavResultsDir, 'JNDStatTable.csv');
writetable(JNDStatTable, JNDTableCSV);
end

function sortStr = initSortedStruct(numCond, numRun)
% sortStr = initSortedStruct(numCond) initializes the structure that will
% store the pooled results for each subject, or group of subjects. It is
% created to have different sizes, based on the number of conditions that
% are being tested against. I think this should generalize to subconditions
% of conditions, or two condition crossing, but I have not tested that, and
% currently the above scripts only consider one condition to test against. 

% Basic info about the session, the recordings, the subjects
sortStr.expType = [];
sortStr.subject = [];
sortStr.gender  = [];
sortStr.age     = [];
sortStr.f0      = [];
sortStr.curSess = [];
sortStr.studyID = [];
sortStr.runs    = [];

sortStr.instructions = [];
sortStr.selectOpt    = [];

sortStr.reversalsReached = [];
sortStr.trialsCompleted  = [];
sortStr.timeElapsed      = [];

sortStr.distProgression = cell(numRun, 1);

sortStr.trialsAtReversals = cell(numRun, 1);
sortStr.distAtReversals   = cell(numRun, 1);

sortStr.trialsAtCorrectOpt1   = cell(numRun, 1);
sortStr.distAtCorrectOpt1     = cell(numRun, 1);
sortStr.trialsAtIncorrectOpt1 = cell(numRun, 1);
sortStr.distAtIncorrectOpt1   = cell(numRun, 1);
sortStr.trialsAtCorrectOpt2   = cell(numRun, 1);
sortStr.distAtCorrectOpt2     = cell(numRun, 1);
sortStr.trialsAtIncorrectOpt2 = cell(numRun, 1);
sortStr.distAtIncorrectOpt2   = cell(numRun, 1);

sortStr.JNDScores         = [];
sortStr.LastSetAccuracies = [];
sortStr.catchAccuracies   = [];

sortStr.obvSubj          = {};
sortStr.obvAge           = [];
sortStr.obvGender        = {};
end

function polRes = combineCondTrials(pA, curRes, polRes)

whichCondAr = strcmp(pA.cond, eval(pA.condVar));
wC          = find(whichCondAr == 1);            % Which Condition?

polRes.runs         = cat(1, polRes.runs, {curRes.run});

polRes.instructions = curRes.instructions;
polRes.selectOpt    = curRes.selectOpt;

polRes.reversalsReached = cat(1, polRes.reversalsReached, curRes.reversalsReached);
polRes.trialsCompleted  = cat(1, polRes.trialsCompleted, curRes.trialsCompleted);
polRes.timeElapsed      = cat(1, polRes.timeElapsed, curRes.timeElapsed);

polRes.JNDScores         = cat(1, polRes.JNDScores, curRes.JNDScore);
polRes.LastSetAccuracies = cat(1, polRes.LastSetAccuracies, curRes.LastSetAccuracy);
polRes.catchAccuracies   = cat(1, polRes.catchAccuracies, curRes.catchAccuracy);

polRes.obvSubj         = cat(1, polRes.obvSubj, curRes.participant);
% polRes.obvAge          = cat(1, polRes.obvAge, curRes.age);
polRes.obvGender       = cat(1, polRes.obvGender, curRes.gender);
end

function sortStruc = meanCondTrials(sortStruc)
JNDScores  = sortStruc.JNDScores;
Accuracies = sortStruc.LastSetAccuracies;

numJNDScores   = length(JNDScores);

sortStruc.JNDScoreMean        = round(mean(JNDScores), 2);
sortStruc.JNDScoreSE          = std(JNDScores)/sqrt(numJNDScores);
sortStruc.lastSetAccuracyMean = round(mean(Accuracies), 1);
sortStruc.lastSetAccuracySE   = std(Accuracies)/sqrt(numJNDScores);
end

function pA = identifyPooledJNDResultsSet()

end