function DRF_SomGenConfig()
% Run this script to generate Pooled Analysis configuration files for 
% Dissertation Project Experiment 1: Somatosensory Feedback Perturbation

project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control'; 
pAnalysis     = 'DRF_Som';

dirs               = dfDirs(project);
dirs.SavResultsDir = fullfile(dirs.Results, 'Pooled Analyses', pAnalysis);
dirs.SavConfigFile = fullfile(dirs.SavResultsDir, [pAnalysis 'PooledConfig.mat']);

if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

% Edit the parts between the lines to modify the pooled analysis variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Participants involved in analysis
pooledParticipants = {'DRF2',...
                      'DRF4',...
                      'DRF5',...
                      'DRF6',...
                      'DRF7',...
                      'DRF9',...
                      'DRF10',...
                      'DRF12',...
                      'DRF13',...
                      'DRF14',...
                      'DRF15',...
                      'DRF16',...
                      'DRF17',...
                      'DRF18',...
                      'DRF19',...
                      'DRF20'};

% Runs for each participant. There should be an equal number of runs for
% each participant. Each row will be the runs to 
pooledRuns  = {'SF1', 'SF2', 'SF3', 'SF4'};
           
% Conditions to test against
testingConditions = {'Voice Feedback', 'Masking Noise'};
pubConditions     = {'notMasked', 'Masked'};

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
cF.pubCond      = pubConditions;

save(dirs.SavConfigFile, 'cF');
fprintf('%s Pooled Analysis Config File Generated!\n', pAnalysis)
end