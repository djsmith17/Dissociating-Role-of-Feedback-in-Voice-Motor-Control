
project = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';

exp = 'NewSensorTest_Winter18';
runs = {'DS5', 'DS6', 'DS7', 'DS8', 'DS9', 'DS10', 'DS11', 'DS12', 'DS13', 'DS14', 'DS15', 'DS16', 'DS17'};
numRuns = length(runs);

dirs = dfDirs(project);

dirs.SavedResultsDir = fullfile(dirs.Results, exp);

f0b = 100;
pF  = 1;
iRF = 0;

for ii = 1:numRuns
    curRun = runs{ii};
    dirs.SavedFileDir = fullfile(dirs.SavData, exp, curRun);
    theFile           = fullfile(dirs.SavedFileDir, [exp curRun 'NSD.mat']);
    
    load(theFile) % returns NSD    

    [niAn, niRes] = dfAnalysisNIDAQ(dirs, NSD.expParam, NSD.DAQin, f0b, 0, iRF, pF);        
    drawDAQAlignedPressure(niRes, dirs.SavResultsDir, 1)
    close all
end