function [auAn, auRes] = dfAnalysisAudapter(dirs, expParam, rawData, bTf0b, AudFlag)
%Analyses the microphone data from the somatosensory perturbation
%experiment. Measures the change in f0 over each trial, and each run for a
%given participant. At the end it approximates a general response to
%inflation to be used in the auditory perturbation experiment

%Requires the Signal Processing Toolbox

%Identify some starting variables
auAn.subject   = expParam.subject;
auAn.run       = expParam.run;
auAn.curSess   = expParam.curSess;
auAn.gender    = expParam.gender;
auAn.AudFB     = expParam.AudFB;
auAn.AudFBSw   = expParam.AudFBSw;
auAn.bTf0b     = bTf0b;
auAn.trialType = expParam.trialType;

fprintf('\nStarting Audapter Analysis for %s, %s\n', auAn.subject, auAn.run)

auAn.sRate    = expParam.sRateAnal;
auAn.trialLen = expParam.trialLen;
auAn.numSamp  = auAn.sRate*auAn.trialLen;
auAn.numTrial = expParam.numTrial;
auAn.expTrigs = expParam.trigs(:,:,1); %Time
auAn.anaTrigs = expParam.trigs(:,:,3);

auAn.time     = (0:1/auAn.sRate:(auAn.numSamp-1)/auAn.sRate)';
auAn.audioM   = [];
auAn.audioH   = [];

dnSamp = 20;
auAn.sRateDN = auAn.sRate/dnSamp;
auAn.timeDN = dnSampleSignal(auAn.time, dnSamp);
auAn.anaTrigsDN = auAn.expTrigs*auAn.sRateDN;

for ii = 1:auAn.numTrial
    data = rawData(ii);
    
    Mraw = data.signalIn;     % Microphone
    Hraw = data.signalOut;    % Headphones
   
    auAn.audProcDel = data.params.frameLen*12;
    [mic, head, saveT, saveTmsg] = preProc(Mraw, Hraw, auAn.sRate, auAn.audProcDel, auAn.expTrigs(ii,1));

%     OST  = data.ost_stat;
    if saveT == 0 %Don't save the trial :(
        fprintf('%s Trial %d not saved. %s\n', auAn.curSess, ii, saveTmsg)
    elseif saveT == 1 %Save the Trial

    end

    auAn.audioM = cat(2, auAn.audioM, mic);
    auAn.audioH = cat(2, auAn.audioH, head);
end

[auAn.ContTrials, auAn.contIdx] = find(auAn.trialType == 0);
[auAn.PertTrials, auAn.pertIdx] = find(auAn.trialType == 1);
auAn.numContTrials = sum(auAn.ContTrials);
auAn.numPertTrials = sum(auAn.PertTrials);

auAn.contTrig = auAn.expTrigs(auAn.contIdx,:);
auAn.pertTrig = auAn.expTrigs(auAn.pertIdx,:);

%The Audio Analysis
f0Flag = 1;
auAn = dfAnalysisAudio(dirs, auAn, AudFlag, f0Flag);

lims  = identifyLimits(auAn);
auRes = packResults(auAn, lims);
end

function [micP, headP, saveT, saveTmsg] = preProc(micR, headR, fs, audProcDel, spanSt)
%This function performs pre-processing on the recorded audio data before
%frequency analysis is applied. This function takes the following inputs:

%micR:       Raw Microphone signal
%headR:      Raw Headphone signal
%fs:         Recording sampling rate
%audProcDel: The delay that results from Audapter processing audio

%This function outputs the following
%micP:     Processed Microphone signal
%headP:    Processed Headphone signal
%saveT:    Boolean toggle to determine if the trial should be saved
%saveTmsg: Reason, if any, that the trial was thrown out

Mfull = micR(1:(end-audProcDel));
Hfull = headR((audProcDel+1):end); 

x = double(Mfull); 
y = double(Hfull);

pp.lenSig = length(x);
pp.t = 0:1/fs:(pp.lenSig-1)/fs;

pp.thresh = 0.3;
[B,A] = butter(4,40/(fs/2)); %Low-pass filter under 40

%Envelope the signal removing all high frequncies. 
%This shows the general change in amplitude over time. 
pp.xenv  = filter(B,A,abs(x));  

%The largest peak in the envelope theoretically occurs during voicing
pp.maxPeak = max(pp.xenv); 

%I don't want to start my signal at the max value, so start lower down on
%the envelope as a threshold
pp.threshIdx     = find(pp.xenv > pp.thresh*pp.maxPeak); 

%The first index of the theoretical useable signal (Voice onset)
pp.voiceOnsetInd = pp.threshIdx(1);  

%The rest of the signal base the first index...are there any dead zones??
pp.fallOffLog = pp.xenv(pp.voiceOnsetInd:end) < pp.thresh*pp.maxPeak;
pp.chk4Break = sum(pp.fallOffLog) > 0.3*fs; %Last longer than 300ms

[B,A]    = butter(4,(300)/(fs/2));
filtx    = filtfilt(B,A,x); %Low-pass filtered under 500Hz
filty    = filtfilt(B,A,y); %Low-pass filtered under 500Hz
 
micP     = filtx; %Take the whole signal for now
headP    = filty; %Same indices as for mic 

micP     = Mfull;
headP    = Hfull;

if pp.t(pp.voiceOnsetInd) > spanSt
    saveT = 0;  
    saveTmsg = 'Participant started too late!!';
elseif pp.chk4Break
    saveT = 0;
    saveTmsg = 'Participant had a voice break!!';
else
    saveT = 1;
    saveTmsg = 'Everything is good'; 
end

pp.saveT    = saveT;
pp.SaveTmsg = saveTmsg;
end

function sensorDN = dnSampleSignal(sensor, dnSamp)
[numSamp, numTrial] = size(sensor);
numSampDN = numSamp/dnSamp;

sensorDN = zeros(numSampDN, numTrial);
for i = 1:numSampDN
    sensorDN(i,:) = mean(sensor((1:dnSamp) + dnSamp*(i-1),:));
end
end

function lims = identifyLimits(An)

%%%%%%%%%%%lims.audio%%%%%%%%%%%
%Full Individual Trials: f0 Audio
if ~isempty(An.audioMf0_pPP)
    pertTrialsM = An.audioMf0_pPP;
    pertTrialsH = An.audioHf0_pPP;
    sec = 100:700;

    uLMa = max(pertTrialsM(sec,:));
    uLMa(find(uLMa > 150)) = 0;
    lLMa = min(pertTrialsM(sec,:));
    lLMa(find(lLMa < -150)) = 0;

    uLM = round(max(uLMa)) + 20;
    lLM = round(min(lLMa)) - 20;
       
    uLHa = max(pertTrialsH(sec,:));
    uLHa(find(uLHa > 150)) = 0;
    lLHa = min(pertTrialsH(sec,:));
    lLHa(find(lLHa < -150)) = 0;
    
    uLH = round(max(uLHa)) + 20;
    lLH = round(min(lLHa)) - 20;
    
    if uLH > uLM
        uLMH = uLH;
    else
        uLMH = uLM;
    end

    if lLH < lLM
        lLMH = lLH;
    else
        lLMH = lLM;
    end
    
    lims.audioM         = [0 4 lLM uLM];
    lims.audioAudRespMH = [0 4 lLMH uLMH];
else
    lims.audioM         = [0 4 -20 20];
    lims.audioAudRespMH = [0 4 -100 100];
end

%%%%%%%%%%%lims.audioMean%%%%%%%%%%%
%Section Mean Pertrubed Trials: f0 Audio 
if ~isempty(An.audioMf0_meanp)
    [~, Imax] = max(An.audioMf0_meanp(:,1)); %Max Pert Onset
    upBoundOn = round(An.audioMf0_meanp(Imax,1) + An.audioMf0_meanp(Imax,2) + 10);
    [~, Imin] = min(An.audioMf0_meanp(:,1)); %Min Pert Onset
    lwBoundOn = round(An.audioMf0_meanp(Imin,1) - An.audioMf0_meanp(Imin,2) - 10);

    [~, Imax] = max(An.audioMf0_meanp(:,3)); %Max Pert Offset
    upBoundOf = round(An.audioMf0_meanp(Imax,3) + An.audioMf0_meanp(Imax,4) + 10);
    [~, Imin] = min(An.audioMf0_meanp(:,3)); %Min Pert Offset
    lwBoundOf = round(An.audioMf0_meanp(Imin,3) - An.audioMf0_meanp(Imin,4) - 10);

    if upBoundOn > upBoundOf
        upBoundM = upBoundOn;
    else
        upBoundM = upBoundOf;
    end

    if lwBoundOn < lwBoundOf
        lwBoundM = lwBoundOn;
    else
        lwBoundM = lwBoundOf;
    end

    lims.audioMean = [-0.5 1.0 lwBoundM upBoundM];
else
    lims.audioMean = [-0.5 1.0 -50 50];
end

%%%%%%%%%%%lims.audioMH%%%%%%%%%%%%
%Section Mean Pertrubed Trials: f0 Audio 
if ~isempty(An.audioHf0_meanp)
    [~, Imax] = max(An.audioHf0_meanp(:,1)); %Max Pert Onset
    upBoundOn = round(An.audioHf0_meanp(Imax,1) + An.audioHf0_meanp(Imax,2) + 10);
    [~, Imin] = min(An.audioHf0_meanp(:,1)); %Min Pert Onset
    lwBoundOn = round(An.audioHf0_meanp(Imin,1) - An.audioHf0_meanp(Imin,2) - 10);

    [~, Imax] = max(An.audioHf0_meanp(:,3)); %Max Pert Offset
    upBoundOf = round(An.audioHf0_meanp(Imax,3) + An.audioHf0_meanp(Imax,4) + 10);
    [~, Imin] = min(An.audioHf0_meanp(:,3)); %Min Pert Offset
    lwBoundOf = round(An.audioHf0_meanp(Imin,3) - An.audioHf0_meanp(Imin,4) - 10);

    if upBoundOn > upBoundOf
        upBoundH = upBoundOn;
    else
        upBoundH = upBoundOf;
    end

    if lwBoundOn < lwBoundOf
        lwBoundH = lwBoundOn;
    else
        lwBoundH = lwBoundOf;
    end
    
    if upBoundH > upBoundM
        upBoundMH = upBoundH;
    else
        upBoundMH = upBoundM;
    end
    
    if lwBoundH < lwBoundM
        lwBoundMH = lwBoundH;
    else
        lwBoundMH = lwBoundM;
    end

    lims.audioMH = [-0.5 1.0 lwBoundMH upBoundMH];
else
    lims.audioMH = [-0.5 1.0 -100 50];
end

end

function res = packResults(auAn, lims)

res.subject = auAn.subject;
res.run     = auAn.run;
res.curSess = auAn.curSess;
res.AudFB   = auAn.AudFB;

res.numTrials     = auAn.numTrial;
res.numContTrials = auAn.numContTrials;
res.numPertTrials = auAn.numPertTrials;
res.contIdx       = auAn.contIdx;
res.pertIdx       = auAn.pertIdx;
res.pertTrig      = auAn.pertTrig;

res.timeA     = auAn.time_audio;
res.f0b       = auAn.f0b;

res.numContTrialsPP = auAn.numContTrialsPP;
res.numPertTrialsPP = auAn.numPertTrialsPP;
res.pertTrigPP      = auAn.pertTrigPP;

%Full Individual Trials: Mic/Head f0 Trace 
res.audioMf0TrialPert = auAn.audioMf0_pPP;
res.audioMf0TrialCont = auAn.audioMf0_cPP;
res.audioHf0TrialPert = auAn.audioHf0_pPP;
res.audioHf0TrialCont = auAn.audioHf0_cPP;
res.limitsA           = lims.audioM;
res.limitsAudRes      = lims.audioAudRespMH;

%Sections Trials: Mic/Head f0
res.secTime          = auAn.secTime;
res.audioMf0SecPert  = auAn.audioMf0_Secp;
res.audioMf0SecCont  = auAn.audioMf0_Secc;
res.audioHf0SecPert  = auAn.audioHf0_Secp;
res.audioHf0SecCont  = auAn.audioHf0_Secc;

%Mean Sectioned Trials: Mic/Head f0 Trace 
res.audioMf0MeanPert = auAn.audioMf0_meanp; % [MeanSigOn 90%CI MeanSigOff 90%CI]
res.audioMf0MeanCont = auAn.audioMf0_meanc;
res.audioHf0MeanPert = auAn.audioHf0_meanp;
res.audioHf0MeanCont = auAn.audioHf0_meanc;
res.limitsAmean      = lims.audioMean;
res.limitsAMH        = lims.audioMH;

%Inflation Response
res.respVar   = auAn.respVar;
res.respVarM  = auAn.respVarMean;
res.respVarSD = auAn.respVarSD;
res.InflaStimVar = auAn.InflaStimVar;
end