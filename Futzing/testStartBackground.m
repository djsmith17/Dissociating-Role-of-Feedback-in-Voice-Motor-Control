function testStartBackground()

close all;
expParam.expType  = 'Somatosensory Perturbation_Perceptual';
expParam.trialLen = 4;
expParam.numTrial = 5;
expParam.perCatch = 1;
expParam.sRateQ   = 8000;
expParam.sRateAnal= 16000;

[s, niCh, nVS]  = initNIDAQ(expParam.trialLen, 'Dev2');

expParam.trialType = dfSetTrialOrder(expParam.numTrial, expParam.perCatch);

[expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType, 1);

expParam.resPause = 3;

% foregroundVS(expParam, s, niCh, nVS);
backgroundVS(expParam, s, niCh, nVS)
end

function foregroundVS(expParam, s, niCh, nVS)

DAQin = []; DAQtime = [];
for ii = 1:expParam.numTrial
    NIDAQsig = [expParam.sigs(:,ii) nVS];
    queueOutputData(s, NIDAQsig);
    fprintf('Running Trial %d\n', ii)
    
    tic
    [data_DAQ, time] = s.startForeground;
    toc

    DAQin   = cat(3, DAQin, data_DAQ);
    DAQtime = cat(3, DAQtime, time);

    pause(expParam.resPause)      
end

figure
plot(DAQtime(:,1,1), DAQin(:,1,1))
end

function backgroundVS(expParam, s, niCh, nVS)

s.IsContinuous = true;

allData = [];
ALLSIGS = []; ALLnVS = [];
for ii = 1:expParam.numTrial
    ALLSIGS = cat(1,ALLSIGS, expParam.sigs(:,ii));
    ALLnVS  = cat(1,ALLnVS, nVS);
end
NIDAQsig = [ALLSIGS ALLnVS];
sd = saveNIDAQ;
    
lrec  = addlistener(s,'DataAvailable', ...
                  @(src,event) updateNIDAQdata(sd, event.TimeStamps, event.Data));

lsend = addlistener(s,'DataRequired', ...
                  @(src,event) src.queueOutputData(NIDAQsig));

queueOutputData(s, NIDAQsig);
% fprintf('Running Trial %d\n', ii)


startBackground(s);
tic
pause(20);
toc
stop(s)
delete(lrec)
delete(lsend)

time = sd.Time;
data = sd.Data(:,1);

figure
plot(time, data)
end

function plotData(src,event)
     plot(event.TimeStamps(:,1), event.Data(:,1))
%      axis([0 4 0 5])
     hold on
end