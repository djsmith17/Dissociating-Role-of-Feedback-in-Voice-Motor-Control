function [auAn, auRes] = dfAnalysisAudapter(dirs, expParam, rawData, DAQin, bTf0b, AudFlag)

%dirs, expParam, DAQin, bTf0b, AudFlag

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

auAn.dnSamp   = 10;
auAn.sRate    = expParam.sRateAnal;
auAn.trialLen = expParam.trialLen;
auAn.numSamp  = expParam.sRate*expParam.trialLen;
auAn.numTrial = expParam.numTrial;
auAn.expTrigs = expParam.trigs(:,:,1); %Time
auAn.anaTrigs = expParam.trigs(:,:,3);

[auAn.ContTrials, auAn.contIdx] = find(auAn.trialType == 0);
[auAn.PertTrials, auAn.pertIdx] = find(auAn.trialType == 1);
auAn.numContTrials = sum(auAn.ContTrials);
auAn.numPertTrials = sum(auAn.PertTrials);

auAn.contTrig = auAn.anaTrigs(:, auAn.contIdx);
auAn.pertTrig = auAn.anaTrigs(:, auAn.pertIdx);

auAn.time     = (0:1/auAn.sRate:(auAn.numSamp-1)/auAn.sRate)';

auAn.actualRecLen = length(rawData(1).signalIn)/auAn.sRate;
auAn.frameT       = linspace(0,auAn.actualRecLen, 2053);
    
auAn.anaInds(:,1) = auAn.winSts;                      %Start indice for analysis based on EvalStep
auAn.anaInds(:,2) = auAn.winSts + auAn.winLenP - 1;   %Stop indice for analysis based on EvalStep
auAn.time         = mean(auAn.anaInds, 2)/auAn.sRate;  %Vector of time points roughly centered on start and stop points of analysis
auAn.baseTimeInd  = auAn.time > 0.5 & auAn.time < 1.0; %The period 0.5s before perturbations

%Analysis Time Steps for Sectioned trial
auAn.winStsSec  = 1:auAn.tStepP:(auAn.totEveLenP-auAn.winLenP); %Starting indices for each analysis window
auAn.numWinSec  = length(auAn.winStsSec); %Number of analysis windows;       

auAn.anaIndsSec(:,1) = auAn.winStsSec;                      % Start indice for analysis based on EvalStep
auAn.anaIndsSec(:,2) = auAn.winStsSec + auAn.winLenP - 1;   % Stop indice for analysis based on EvalStep
auAn.timeSec         = mean(auAn.anaIndsSec, 2)/auAn.sRate; % Vector of time points roughly centered on start and stop points of analysis

auRes.time          = auAn.time;
auRes.timeSec       = auAn.timeSec;
auRes.runTrialOrder = [];
auRes.allTrialf0    = [];
auRes.allTrialf0_St = [];
auRes.allTrialf0_Sp = [];
auRes.allTrialf0b   = [];
auRes.allTrialForce = [];
auRes.allTrialTrigs = auAn.trigsT;
auRes.trialCount    = [];
auRes.meanTrialf0_St = [];
auRes.meanTrialf0_Sp = [];
auRes.meanTrialf0b   = [];
auRes.meanTrialForce_St = [];
auRes.meanTrialForce_Sp = [];

for ii = 1:auAn.numTrial
    data = rawData(ii);
    
    Mraw = data.signalIn;     % Microphone
    Hraw = data.signalOut;    % Headphones
    OST  = data.ost_stat;
    audProcDel = data.params.frameLen*4;
    
    [mic, head, saveT, saveTmsg] = preProc(Mraw, Hraw, auAn.sRate, audProcDel, auAn.trigsT(ii,1));

    if saveT == 0 %Don't save the trial :(
        fprintf('%s Trial %d not saved. %s\n', auAn.curSess, ii, saveTmsg)
    elseif saveT == 1 %Save the Trial
%         fprintf('%s Trial %d saved\n', auAn.curSess, ii)
        
        %Pert Full Trial
        Trialf0Raw = signalFrequencyAnalysis(mic, head, auAn.sRate, auAn, [], 0);
              
        %Pert Section: Start
        Trialf0Raw_St = signalFrequencyAnalysis(mic, head, auAn.sRate, auAn, auAn.trigsA(ii,1), 1);
        %Pert Section: Stop
        Trialf0Raw_Sp = signalFrequencyAnalysis(mic, head, auAn.sRate, auAn, auAn.trigsA(ii,2), 1);
        
%         trig = auAn.trigsT(ii,1);
%         trigSt = trig - 0.5;
%         trigSp = trig + 1.0; 
%         tFramesSt = find(trigSt >= auAn.frameT); tFrameSt = tFramesSt(end);
%         tFramesSp = find(trigSp >= auAn.frameT); tFrameSp = tFramesSp(end);
%         frameT_Start = auAn.frameT(tFrameSt:tFrameSp);
%         OST_Start = OST(tFrameSt:tFrameSp);
%         headedup = find(OST_Start == 3); headUp = headedup(1);
%         whenitactuallystarted = frameT_Start(headUp);
%         whatTHEDIFF = whenitactuallystarted - trig
        

%         prePertInd = auAn.time < 0.5;                    % Grab the first 0.5s, should be no stimulus
        f0b = round(mean(Trialf0Raw(auAn.baseTimeInd, 1))); % Baseline fundamental frequency of mic data
        
        auRes.runTrialOrder = cat(1, auRes.runTrialOrder, auAn.trialType(ii));
        auRes.allTrialf0 = cat(3, auRes.allTrialf0, Trialf0Norm);
        auRes.allTrialf0_St = cat(3, auRes.allTrialf0_St, Trialf0Norm_St);
        auRes.allTrialf0_Sp = cat(3, auRes.allTrialf0_Sp, Trialf0Norm_Sp);
        auRes.allTrialf0b   = cat(1, auRes.allTrialf0b, f0b);            %Baseline fundamental frequencies
        auRes.allTrialForce = cat(3, auRes.allTrialForce, TrialForce);   %Force sensor values;
    end
end

%Sort trials within a given run by trial type and find average across trials
[auRes.meanTrialf0_St, auRes.meanTrialForce_St, auRes.trialCount] = sortTrials(auRes.allTrialf0_St, auRes.allTrialForce, auRes.runTrialOrder);
[auRes.meanTrialf0_Sp, auRes.meanTrialForce_Sp, auRes.trialCount] = sortTrials(auRes.allTrialf0_Sp, auRes.allTrialForce, auRes.runTrialOrder);
auRes.meanTrialf0b = round(mean(auRes.allTrialf0b,1));

%The Audio Analysis
auAn = dfAnalysisAudio(dirs, auAn, AudFlag);

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

mic  = micR(1:(end-audProcDel));
head = headR((audProcDel+1):end); 

x = double(mic); 
y = double(head);

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

function Trialf0ResultsRaw = signalFrequencyAnalysis(mic, head, fs, auAn, trig, sec)
%Finds the change in fundamental frequency of windowed signal

%Inputs
%mic:  post-processed single-trial Microphone signal
%head: post-processed signle-trial Headphone signal
%fs:   sampling frequency of audio signals
%auAn: audapter analysis structure containing other relevant var
%trig: trigger points in audio signals when event occurs (onset/offset)
%sec:  boolean to check for full signal analysis or sectioned analysis

%Outputs:
%Trialf0ResultsRaw: Array with rows equal to the number of windows. 
%                   The first column is the time value that the window is 
%                   centered around. 
%                   The second column is the fundamental frequency of the 
%                   windowed microphone signal. 
%                   The third column is the fundamental frequency of the 
%                   windowed headphone signal.

St = trig - auAn.preEveLenP; 
Sp = trig + auAn.posEveLenP - 1;

%Grab a big chuck of the signal centered around the event
if sec == 1
    try
        mic = mic(St:Sp);
        head = head(St:Sp);
    catch
        disp('Sp was too long yo!')
        numSamp = length(mic);
        mic = mic(St:numSamp);
        head = head(St:numSamp);
    end
    
    numWin  = auAn.numWinSec;
    anaInds = auAn.anaIndsSec; 
else
    numWin  = auAn.numWin;
    anaInds = auAn.anaInds;    
end

Trialf0ResultsRaw = [];
for ii = 1:numWin
    startPt  = anaInds(ii,1);
    stopPt   = anaInds(ii,2);

    mic_win   = mic(startPt:stopPt);
    head_win  = head(startPt:stopPt);
    
    f0_M = dfCalcf0Chile(mic_win,fs);
    f0_H = dfCalcf0Chile(head_win,fs);
    
%     [f0_time, f0_value, SHR, f0_candidates] = shrp(mic_win, fs);
    if f0_M < 50 || f0_M > 300
        fprintf('I calculated a f0 of %d. Replacing it.\n', f0_M)
        f0_M = 0; %Trialf0ResultsRaw(ii-1,2);
    end
    
    if f0_H < 50 || f0_H > 300
        disp('I had some difficulty calculating f0_H')
        f0_H = 0;
    end
    
    Trialf0ResultsRaw = cat(1, Trialf0ResultsRaw, [f0_M f0_H]);
end
end

function [meanTrialf0, meanTrialForce, trialCount] = sortTrials(allTrialf0, allTrialForce, runTrialOrder)
%This function separates the trials by control or catch trials and
%finds the mean f0 trace and 95% Confidence Interval over multiple trials 
%of a type. 

%Sort Trials by type 
PertVals  = unique(runTrialOrder(:,1));
nPertVals = length(PertVals);

meanTrialf0    = [];
meanTrialForce = [];
trialCount   = zeros(1, nPertVals);
for i = 1:nPertVals
    ind   = runTrialOrder == PertVals(i);
    nType = sum(ind);
    
    TrialsofaType_f0  = allTrialf0(:,:,ind);    
    micMean_f0   = mean(squeeze(TrialsofaType_f0(:,1,:)),2);
    headMean_f0  = mean(squeeze(TrialsofaType_f0(:,2,:)),2);
    
    micSTD_f0    = std(squeeze(TrialsofaType_f0(:,1,:)),0,2);
    headSTD_f0   = std(squeeze(TrialsofaType_f0(:,2,:)),0,2);    
    
    SEM_f0 = micSTD_f0/sqrt(nType); % Standard Error
    CIM_f0 = 1.96*SEM_f0; % 95% confidence Interval
    
    SEH_f0 = headSTD_f0/sqrt(nType); % Standard Error
    CIH_f0 = 1.96*SEH_f0; % 95% confidence Interval
                    
%     ts  = tinv([0.025  0.975],numT-1);      % T-Score
%     CIM = micMean + ts*SEM; 
               
%     ts  = tinv([0.025  0.975],numT-1);      % T-Score
%     CIH = headMean + ts*SEM; 
   f0resultSet  = [micMean_f0 CIM_f0 headMean_f0 CIH_f0];
   meanTrialf0  = cat(3, meanTrialf0, f0resultSet);
   
   TrialsofaType_Force  = allTrialForce(:,:,ind);
   collMean_F = mean(squeeze(TrialsofaType_Force(:,1,:)),2);
   neckMean_F = mean(squeeze(TrialsofaType_Force(:,2,:)),2);
   
   collSTD_F  = std(squeeze(TrialsofaType_Force(:,1,:)),0,2);
   neckSTD_F  = std(squeeze(TrialsofaType_Force(:,2,:)),0,2);
   
   SEC_F = collSTD_F/sqrt(nType); % Standard Error
   CIC_F = 1.96*SEC_F; % 95% confidence Interval
    
   SEN_F = neckSTD_F/sqrt(nType); % Standard Error
   CIN_F = 1.96*SEN_F; % 95% confidence Interval 
   
   ForceresultSet = [collMean_F CIC_F neckMean_F CIN_F];
   meanTrialForce = cat(3, meanTrialForce, ForceresultSet);
   
   trialCount(i)= nType;
end
end

function lims = identifyLimits(An)

%%%%%%%%%%%lims.audio%%%%%%%%%%%
%Full Individual Trials: f0 Audio
if ~isempty(An.audioMf0_pPP)
    pertTrials = An.audioMf0_pPP;
    sec = 100:700;

    alluL = max(pertTrials(sec,:));
    alluL(find(alluL > 150)) = 0;
    alllL = min(pertTrials(sec,:));
    alllL(find(alllL < -150)) = 0;

    uL = round(max(alluL)) + 20;
    lL = round(min(alllL)) - 20;
    lims.audio      = [0 4 lL uL];
else
    lims.audio      = [0 4 -20 20];
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
        upBoundSec = upBoundOn;
    else
        upBoundSec = upBoundOf;
    end

    if lwBoundOn < lwBoundOf
        lwBoundSec = lwBoundOn;
    else
        lwBoundSec = lwBoundOf;
    end

    lims.audioMean = [-0.5 1.0 lwBoundSec upBoundSec];
else
    lims.audioMean = [-0.5 1.0 -50 50];
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

res.timeS      = auAn.time_DN;
res.sensorP    = auAn.sensorP_p; %Individual Processed perturbed trials. 
res.lagTimeP   = auAn.lagsPres;
res.lagTimePm  = auAn.meanLagTimeP;
res.riseTimeP  = auAn.riseTimeP;
res.riseTimePm = auAn.riseTimePm;
res.OnOfValP   = auAn.OnOfValP;
res.OnOfValPm  = auAn.OnOfValPm;
res.limitsP    = lims.pressure;

res.timeSAl   = auAn.time_Al;
res.sensorPAl = auAn.sensorP_Al;
res.limitsPAl = lims.pressureAl;

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
res.limitsA           = lims.audio;

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

%Inflation Response
res.respVar   = auAn.respVar;
res.respVarM  = auAn.respVarMean;
res.respVarSD = auAn.respVarSD;
res.InflaStimVar = auAn.InflaStimVar;
end