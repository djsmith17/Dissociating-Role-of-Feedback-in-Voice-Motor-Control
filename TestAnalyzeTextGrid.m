clear all;
FileDir = 'C:\users\djsmith\desktop';
FileNm  = 'Pilot22AF1Trial4_micIn.TextGrid';

File = fullfile(FileDir, FileNm);

[tg, fid] = tgRead(File);

nTiers = tgGetNumberOfTiers(tg);
tierIndex = tgI(tg, 'words');

numIntervals = tgGetNumberOfIntervals(tg, tierIndex);

for ii = 1:numIntervals
    dur = tgGetIntervalDuration(tg, tierIndex, ii)
    name = tgGetLabel(tg, tierIndex, ii)
end

% tgIsPointTier(tg, 1)

figure
tgPlot(tg)