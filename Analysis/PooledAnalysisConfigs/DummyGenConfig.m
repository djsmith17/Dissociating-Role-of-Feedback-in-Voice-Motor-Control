function DummyGenConfig()
% Run this script to generate Pooled Analysis configuration files for 
% Masking Noise piloting

project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control'; 
pAnalysis     = 'Dummy';

dirs               = dfDirs(project);
dirs.SavResultsDir = fullfile(dirs.Results, 'Pooled Analyses', pAnalysis);
dirs.SavConfigFile = fullfile(dirs.SavResultsDir, [pAnalysis 'PooledConfig.mat']);

if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavRdesultsDir)
end

% Edit the parts between the lines to modify the pooled analysis variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Participants involved in analysis
pooledParticipants = {'DRF_MN18'};

% Runs for each participant. There should be an equal number of runs for
% each participant. Each row will be the runs to 
pooledRuns  = {'SF1', 'SF2', 'SF3'};
           
% Conditions to test against
testingConditions = {'Voice Feedback', 'AC Masking Noise', 'AC/BC Masking Noise'};

% The Recording Variable to check for the condition
condVar = 'curRes.AudFB'; 

% How do you want to title the Result Plots?
testExt    = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cF.participants = pooledParticipants;
cF.runs         = pooledRuns;
cF.cond         = testingConditions;
cF.condVar      = condVar;
cF.testExt      = testExt;

save(dirs.SavConfigFile, 'cF');
fprintf('%s Pooled Analysis Config File Generated!\n', pAnalysis)
end