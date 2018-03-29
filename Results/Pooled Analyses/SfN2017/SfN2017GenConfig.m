function SfN2017GenConfig()
% Run this script to generate Pooled Analysis configuration files for 
% SfN2017

project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control'; 
pAnalysis     = 'SfN2017';

dirs               = dfDirs(project);
dirs.SavResultsDir = fullfile(dirs.Results, 'Pooled Analyses', pAnalysis);
dirs.SavConfigFile = fullfile(dirs.SavResultsDir, [pAnalysis 'PooledConfig.mat']);

if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

% Edit the parts between the lines to modify the pooled analysis variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Participants involved in analysis
pooledParticipants = {'Pilot24';...
                      'Pilot25';...
                      'Pilot26';...
                      'Pilot22'};

% Runs for each participant. There should be an equal number of runs for
% each participant. Each row will be the runs to 
pooledRuns  = {'SF1', 'SF2', 'SF3', 'SF4';...
               'SF1', 'SF2', 'SF3', 'SF4';...
               'SF1', 'SF2', 'SF3', 'SF4';...
               'SF1', 'SF2', 'SF3', 'SF4'};
           
% Conditions to test against
testingConditions = {'Masking Noise'; 'Normal Voicing'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numPart = length(pooledParticipants);

[p, ~] = size(pooledRuns);

if p ~= numPart
    fprintf('Please double check your Pooled Config File\n')
    fprintf('Number of rows in pooledRuns does not match number of participants\n')
    return
end

cF.participants = pooledParticipants;
cF.runs         = pooledRuns;
cF.cond         = testingConditions;

save(dirs.SavConfigFile, 'cF');
fprintf('%s Pooled Analysis Config File Generated!\n', pAnalysis)
end