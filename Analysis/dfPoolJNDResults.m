function dfPoolJNDResults()

close all
pA.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control'; 
pA.pAnalysis     = 'DRF_JND'; % Change this name to load different pooled data sets Ex: SfN2017, LarynxPos

dirs               = dfDirs(pA.project);
dirs.SavResultsDir = fullfile(dirs.Results, 'Pooled Analyses', pA.pAnalysis);
dirs.PooledConfigF = fullfile(dirs.SavResultsDir, [pA.pAnalysis 'PooledConfig.mat']);

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
        dirs.SavFileDir  = fullfile(dirs.Results, participant, run);                      % Where results are saved
        dirs.SavFile     = fullfile(dirs.SavFileDir, [participant run 'ResultsDRF.mat']); % Run Results file to load

        if exist(dirs.SavFile, 'file') == 0
            error('ERROR: The Results file for Run %s does not exist yet\n', run)
        else   
            load(dirs.SavFile)
            % Returns a results struture of 'res'
        end
        
        subjRes = cat(2, subjRes, res);
    end
    allDataStr = cat(1, allDataStr, subjRes);
end

end

function pA = identifyPooledJNDResultsSet()

end