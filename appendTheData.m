function appendTheData

close all
project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
participant   = 'DRF_ENP1';    % List of multiple participants.
run           = 'SFL1';
trials        = 5;
baselineFile  = 'BV1';            % Baseline Voice information
debug         = 0;

dirs               = dfDirs(project);

dirs.baselineData  = fullfile(dirs.SavData, participant, baselineFile, [participant baselineFile 'DRF.mat']); % Where to find data
      
% Look for the baseline voice info for this participant, then load it
if exist(dirs.baselineData, 'file') == 0
    fprintf('ERROR: Could not find baseline data set at %s\n', dirs.baselineData)
    return
else
    fprintf('Loading baseline data set for %s %s\n', participant, baselineFile)
    load(dirs.baselineData) % Returns DRF
    bV = DRF;
end


for j = 1:trials
    dirs.SavFileDir    = fullfile(dirs.SavData, participant, run, [participant run 'Trial' num2str(j) 'DRF.mat']);  % Where to find data        

    % Look for the recording sessoin raw data for this participant, then load it
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%\n')
    if exist(dirs.SavFileDir, 'file') == 0
        fprintf('ERROR: Could not find saved data set at %s\n', dirs.SavFileDir)
        return
    else
        fprintf('Loading saved data set for %s %s\n', participant, run)
        load(dirs.SavFileDir) % Returns DRF
    end
    
    DRF.expParam.rmsB    = bV.expParam.rmsB;
    DRF.expParam.age     = bV.expParam.age;
    DRF.expParam.gender  = bV.expParam.gender;
    DRF.expParam.f0b     = bV.qRes.meanf0;
    DRF.expParam.targRMS = bV.qRes.meanRMS;
    DRF.expParam.AudFB   = DRF.expParam.AudFB{1};
        
    save(dirs.SavFileDir, 'DRF')
end

end