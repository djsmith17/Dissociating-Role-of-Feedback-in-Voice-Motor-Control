function LarynxPosGenConfig()
% Run this script to generate Pooled Analysis configuration files for 
% Larynx Pos piloting

project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control'; 
pAnalysis     = 'LarynxPos';

dirs               = dfDirs(project);
dirs.SavResultsDir = fullfile(dirs.Results, 'Pooled Analyses', pAnalysis);
dirs.SavConfigFile = fullfile(dirs.SavResultsDir, [pAnalysis 'PooledConfig.mat']);

if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

% Edit the parts between the lines to modify the pooled analysis variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Participants involved in analysis
pooledParticipants = {'Pilot29';...
                      'Pilot30';...
                      'Pilot31';...
                      'Pilot32';...
                      'Pilot33';...
                      'Pilot21';...
                      'Pilot28';...
                      'Pilot0'};

% Runs for each participant. There should be an equal number of runs for
% each participant. Each row will be the runs to 
pooledRuns  = {'DS1', 'DS2', 'DS3', 'DS4', 'DS5', 'DS6'};
           
% Conditions to test against
testingConditions = {'LP', 'uLP', 'CC'};

% The Recording Variable to check for the condition
condVar = 'pA.cond{ceil(str2double(curRes.run(end))/2)}'; 

% How do you want to title the Result Plots?
testExt    = 'Resp_CollarLoc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cF.participants = pooledParticipants;
cF.runs         = pooledRuns;
cF.cond         = testingConditions;
cF.condVar      = condVar;
cF.testExt      = testExt;

save(dirs.SavConfigFile, 'cF');
fprintf('%s Pooled Analysis Config File Generated!\n', pAnalysis)
end