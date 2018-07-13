function rmsB = findRMS2dBRatio(DRF, targdB)

rawData = DRF.rawData;

rms = rawData(1).rms(:,3);
rms = rms(rms ~= 0);

rmsB = rms./(10^(targdB/20));

rmsB = mean(rmsB);
end