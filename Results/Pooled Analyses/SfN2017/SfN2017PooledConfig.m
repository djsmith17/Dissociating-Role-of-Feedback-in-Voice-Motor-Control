
% Participants involved in analysis
pooledParticipants = {'Pilot24'; 'Pilot25'; 'Pilot26'; 'Pilot22'};

% Runs for each participant. There should be an equal number of runs for
% each participant. Each row will be the runs to 
pooledRuns  = {'SF1', 'SF2', 'SF3', 'SF4';...
               'SF1', 'SF2', 'SF3', 'SF4';...
               'SF1', 'SF2', 'SF3', 'SF4';...
               'SF1', 'SF2', 'SF3', 'SF4'};
           
numPart = length(pooledParticipants);

[p, numRuns] = size(pooledRuns);

if p ~= numPart
    fprintf('Please double check your Pooled Config File\n')
    fprintf('Number of rows in pooledRuns does not match number of participants\n')
    return
end