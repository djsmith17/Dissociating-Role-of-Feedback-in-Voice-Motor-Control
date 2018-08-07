function dfInspectRawData(participant, run)

close all
iRD.project     = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
iRD.participant = participant;
iRD.run         = run;

dirs            = dfDirs(iRD.project);
dirs.LoadData   = dirs.RecData;

dirs.SavWavDir  = fullfile(dirs.LoadData, iRD.participant, iRD.run, 'wavFiles'); % Where to find data
    
iRD.micFiles  = dir([dirs.SavWavDir '\*micIn.wav']);
iRD.headFiles = dir([dirs.SavWavDir '\*headOut.wav']);
iRD.numRuns   = length(iRD.micFiles);

if iRD.numRuns == 0
    error('No wav files found for this participant and run')
end

[iW, sP] = setInspectionWindow(iRD);
for ii = 1:iRD.numRuns
   
    curFile = iRD.micFiles(ii).name;
    dirs.SavWavFile = fullfile(dirs.SavWavDir, curFile);
    [trialM, fs] = audioread(dirs.SavWavFile);
    trialLen  = length(trialM);
    trialTime = trialLen/fs;
    time      = linspace(0, trialTime, trialLen);
    
    axes(sP(ii))
    plot(time, trialM)
    title(['Trial ' num2str(ii)])
    axis([0 5 -0.5 0.5]); box off
    if ii ~= 1
        set(gca,'XTickLabel',[],...
                'YTickLabel',[]);
    end     
end
end

function [iW, sP] = setInspectionWindow(iRD)

plotpos = [1680 0];
plotdim = [1680 1050];
iW = figure('Color', [1 1 1]);
set(iW, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

numRow = 4;
numCol = ceil(iRD.numRuns/numRow);

sP = tight_subplot(numRow, numCol, [0.05 0.01],[0.05 0.05],[0.02 0.00]);

suptitle([iRD.participant ' ' iRD.run])
end