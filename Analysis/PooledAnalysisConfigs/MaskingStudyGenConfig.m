function MaskingStudyGenConfig()
% Run this script to generate Pooled Analysis configuration files for 
% Masking Noise piloting

project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control'; 
pAnalysis     = 'MaskingStudy';

dirs               = dfDirs(project);
dirs.SavResultsDir = fullfile(dirs.Results, 'Pooled Analyses', pAnalysis);
dirs.SavConfigFile = fullfile(dirs.SavResultsDir, [pAnalysis 'PooledConfig.mat']);

if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

% Edit the parts between the lines to modify the pooled analysis variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Participants involved in analysis
pooledParticipants = {'DRF_MN2'};

% Runs for each participant. There should be an equal number of runs for
% each participant. Each row will be the runs to 
pooledRuns  = {'SF1', 'SF2', 'SF3'};
           
% Conditions to test against
testingConditions = {'Voice Feedback', 'AC Masking Noise', 'AC/BC Masking Noise'};

% The Recording Variable to check for the condition
condVar = 'curRes.AudFB'; 

% How do you want to title the Result Plots?
pltNameTop = 'MaskingStudy';
testType   = '';
% MaskvVoice Individual
pltNameMVi = {[pltNameTop 'Pilot28' testType]};

pltNameMVm =  [pltNameTop 'MeanSubj' testType];

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
cF.condVar      = condVar;
cF.pltNameMVi   = pltNameMVi;
cF.pltNameMVm   = pltNameMVm;

save(dirs.SavConfigFile, 'cF');
fprintf('%s Pooled Analysis Config File Generated!\n', pAnalysis)
end