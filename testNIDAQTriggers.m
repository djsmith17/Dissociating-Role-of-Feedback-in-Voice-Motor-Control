function testNIDAQTriggers()

close all;
expParam.expType  = 'Somatosensory Perturbation_Perceptual';
expParam.trialLen = 4.1;
expParam.numTrial = 1;
expParam.perCatch = 1;
expParam.sRateQ   = 8000;
expParam.sRateAnal= 16000;
expParam.resPause = 3;

recordedRMS = (5+rand(1991,1))/1000;

targRMS   = 50; %dB just to test
boundsRMS = 3; %+/- dB

expParam.trialType              = dfSetTrialOrder(expParam.numTrial, expParam.perCatch);
[expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType, 1);

close all
[anMsr, H1, H2, H3, fbLines, rec, trigCirc] = dfSetVisFB(targRMS, boundsRMS);

%Initialize NIDAQ
[s, sd, niCh, nVS]  = initNIDAQSess('Dev2', expParam.trialLen);

backgroundVS(expParam, s, niCh, nVS, H2, H3, trigCirc)
end

function backgroundVS(expParam, s, niCh, nVS, H2, H3, trigCirc)

allData = [];
ALLSIGS = []; ALLnVS = [];
for ii = 1:expParam.numTrial
    ALLSIGS = cat(1,ALLSIGS, expParam.sigs(:,ii));
    ALLnVS  = cat(1,ALLnVS, nVS);
end
NIDAQsig = [ALLSIGS ALLnVS];
    
lrec  = addlistener(s,'DataAvailable', ...
                  @(src,event) updateNIDAQdata(sd, event.TimeStamps, event.Data));

lsend = addlistener(s,'DataRequired', ...
                  @(src,event) src.queueOutputData(NIDAQsig)); % @myFunction(src,evt,NIDAQsig)

queueOutputData(s, NIDAQsig);
% fprintf('Running Trial %d\n', ii)

startBackground(s);

pause(4.0)
set(H3,'Visible','off');
set([H2 trigCirc],'Visible','on');
pause(4.0)
set([H2 trigCirc],'Visible','off');

% tic
% pause(20);
% toc
stop(s)
delete(lrec)
delete(lsend)

time = sd.Time;
data = sd.Data(:,1);

close all
figure
plot(time, data)
end