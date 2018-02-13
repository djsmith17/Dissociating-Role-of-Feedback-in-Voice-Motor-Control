function dfQuickAnalysisPlot(DRF)
dirs     = DRF.dirs;
expParam = DRF.expParam;
rawData  = DRF.rawData;

An.subject  = expParam.subject;
An.run      = expParam.run;
An.curSess  = expParam.curSess;
An.f0AnaFile = [An.subject An.run 'f0AnalysisA.mat'];
An.bTf0b    = expParam.bTf0b;
An.sRate    = expParam.sRateAnal;
An.numSamp  = expParam.trialLen*expParam.sRateAnal;
An.numTrial = expParam.numTrial;
An.pertIdx  = find(expParam.trialType == 1);
An.contIdx  = find(expParam.trialType == 0);
An.pertTrig = expParam.trigs(An.pertIdx, :, 1);
An.contTrig = expParam.trigs(An.contIdx, :, 1);

An.time   = (0:1/An.sRate:(An.numSamp-1)/An.sRate)';
An.audioM = []; An.audioH = [];
for i = 1:numTrial 
    mic  = rawData(i).signalIn;
    head = rawData(i).signalOut;
    An.audioM = cat(2, auAn.audioM, mic(An.numSamp));
    An.audioH = cat(2, auAn.audioH, head(An.numSamp));
end

AudFlag = 1;
iRF     = 1; f0Flag  = 1;
An = dfAnalysisAudio(dirs, An, AudFlag, iRF, f0Flag);


drawDAQMeanTrialMicf0(auRes, dirs.SavResultsDir)
end