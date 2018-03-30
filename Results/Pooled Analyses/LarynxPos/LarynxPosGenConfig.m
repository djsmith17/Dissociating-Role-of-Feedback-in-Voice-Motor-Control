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
pooledRuns  = {'DS1', 'DS2', 'DS3', 'DS4', 'DS5', 'DS6';...
               'DS1', 'DS2', 'DS3', 'DS4', 'DS5', 'DS6';...
               'DS1', 'DS2', 'DS3', 'DS4', 'DS5', 'DS6';...
               'DS1', 'DS2', 'DS3', 'DS4_2', 'DS5', 'DS6';...
               'DS1', 'DS2', 'DS3', 'DS4', 'DS5', 'DS6';...
               'DS1', 'DS2', 'DS3', 'DS4', 'DS5', 'DS6';...
               'DS1', 'DS2', 'DS3', 'DS4', 'DS5', 'DS6';...
               'DS1', 'DS2', 'DS3', 'DS4', 'DS5', 'DS6'};
           
% Conditions to test against
testingConditions = {'Masking Noise'; 'Voice Not Shifted'};

% The Recording Variable to check for the condition
condVar = 'res.AudFB'; 

% How do you want to title the Result Plots?
pltNameTop = 'LarynxPos ';
% MaskvVoice Individual
pltNameMVi = {[pltNameTop 'Pilot29' 'Resp_MaskvVoice'],...
              [pltNameTop 'Pilot30' 'Resp_MaskvVoice'],...
              [pltNameTop 'Pilot31' 'Resp_MaskvVoice'],...
              [pltNameTop 'Pilot32' 'Resp_MaskvVoice'],...
              [pltNameTop 'Pilot33' 'Resp_MaskvVoice'],...
              [pltNameTop 'Pilot21' 'Resp_MaskvVoice'],...
              [pltNameTop 'Pilot28' 'Resp_MaskvVoice'],...
              [pltNameTop 'Pilot0' 'Resp_MaskvVoice']};

pltNameMVm =  [pltNameTop 'MeanSubjResp_MaskvVoice'];

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