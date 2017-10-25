function testNIDAQTriggers()

%Uses
%dfSetTrialOrder
%dfMakePertSignal
%initNIDAQSess

close all; clc
expParam.expType  = 'Somatosensory Perturbation_Perceptual';
expParam.dev      = 'Dev2';
expParam.numTrial = 3;
expParam.perCatch = 1;
expParam.trialLen = 4;
expParam.sRateQ   = 8000;
expParam.sRateAnal= 16000;
expParam.resPause = 3;

% recordedRMS = (5+rand(1991,1))/1000;

targRMS   = 50; % dB just to test
boundsRMS = 3;  % +/- dB

expParam.trialType              = dfSetTrialOrder(expParam.numTrial, expParam.perCatch);
[expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType, 1);

[anMsr, H1, H2, H3, fbLines, rec, trigCirc] = dfSetVisFB(targRMS, boundsRMS);

%Initialize NIDAQ
[s, sd, niCh]  = initNIDAQSess(expParam.dev, expParam.numTrial, expParam.trialLen, expParam.sigs);

backgroundVS(expParam, s, sd, H2, H3, trigCirc)
end

function backgroundVS(expParam, s, sd, H2, H3, trigCirc)

queuePertOutputData(sd, s)
startBackground(s);
pause(4.0)
set(H3,'Visible','off');

for ii = 1:expParam.numTrial
    fprintf('Running Trial %d\n', ii)
    set([H2 trigCirc],'Visible','on');
    pause(expParam.trialLen)
    set([H2 trigCirc],'Visible','off');
    pause(expParam.resPause)
end
close all
stop(s)

time = sd.recTime;
data = sd.recData;

figure
plot(time, data)
end