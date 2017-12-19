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
auAn.numSamp  = 
auAn.numTrial = expParam.numTrial;

auAn.trigsT   = expParam.trigs(:,:,1);  %Pregenerated start and stop times for time-alignment with audio data
auAn.trigsA   = expParam.trigs(:,:,3);  %Pregenerated start and stop points (Audapter) for time-alignment with audio data
auAn.trigsQ   = expParam.trigs(:,:,2);  %Pregenerated start and stop points (NIDAQ) for time-alignment with audio data

[auAn.ContTrials, auAn.contIdx] = find(auAn.trialType == 0);
[auAn.PertTrials, auAn.pertIdx] = find(auAn.trialType == 1);
auAn.numContTrials = sum(auAn.ContTrials);
auAn.numPertTrials = sum(auAn.PertTrials);


auAn.winLen  = 0.05;                   % Analysis window length (seconds)
auAn.winLenP = auAn.winLen*auAn.sRate; % Analysis window length (points)
auAn.pOV     = 0.60;                   % Window overlap percentage (decimal)
auAn.tStepP  = auAn.winLenP*(1-auAn.pOV); %Number of points between each analysis window starting indice (Changes with Percent of overlap)
auAn.tStep   = auAn.tStepP/auAn.sRate;    %Amount of time(s) between each analysis window           

auAn.actualRecLen = length(rawData(1).signalIn)/auAn.sRate;
auAn.frameT       = linspace(0,auAn.actualRecLen, 2053);

auAn.recLen  = expParam.trialLen;
auAn.recLenP = expParam.trialLen*auAn.sRate;

auAn.preEveLen  = 0.5; %Amount of time in seconds of observation period before event (onset/offset)
auAn.posEveLen  = 1.0; %Amount of time in seconds of observation period after event (onset/offset)
auAn.totEveLen  = auAn.preEveLen + auAn.posEveLen; %Total length (seconds) of observation time

auAn.preEveLenP = round(auAn.preEveLen*auAn.sRate);  %Amount of points of observation period before event (onset/offset) for Audapter signal
auAn.posEveLenP = round(auAn.posEveLen*auAn.sRate);  %Amount of points of observation period after event (onset/offset) for Audapter signal
auAn.totEveLenP = auAn.preEveLenP + auAn.posEveLenP; %Total length (points) of observation time for Audapter signal
      
auAn.preEveLenQ = round(auAn.preEveLen*auAn.sRateQ);  %Amount of points of observation period before event (onset/offset) for NIDAQ signal
auAn.posEveLenQ = round(auAn.posEveLen*auAn.sRateQ);  %Amount of points of observation period after event (onset/offset) for NIDAQ signal
auAn.totEveLenQ = auAn.preEveLenQ + auAn.posEveLenQ; %Total length (points_NIDAQ) of observation time
auAn.timeQ      = (0:1:(auAn.totEveLenQ-1))/auAn.sRateQ; %Time points_NIDAQ roughly center of start and stop points of analysis

%Analysis Time Steps for Full trial
auAn.winSts  = 1:auAn.tStepP:(auAn.recLenP-auAn.winLenP); %Starting indices for each analysis window
auAn.numWin  = length(auAn.winSts); %Number of analysis windows;       

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

        Trialf0Norm    = normf0(Trialf0Raw, f0b);    % Coverted to cents and normalized        
        Trialf0Norm_St = normf0(Trialf0Raw_St, f0b); % Coverted to cents and normalized
        Trialf0Norm_Sp = normf0(Trialf0Raw_Sp, f0b); % Coverted to cents and normalized
        
        TrialForce = forceSensorAnalysis(DAQin, auAn.trigsQ(ii,1), auAn.sRateQ, auAn); %At the moment only voltage

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

lims = identifyLimits(auRes);

auRes.f0Limits         = [0 auAn.recLen lims.lwBoundSec lims.upBoundSec];
auRes.f0LimitsSec      = [0 auAn.totEveLen lims.lwBoundSec lims.upBoundSec];
auRes.InflaRespLimits  = [0 0.3 -80 10];
auRes.ForceLimits      = [0 auAn.totEveLen 1 3.5];
auRes.PressureLimits   = [0 auAn.totEveLen 20 30];
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

function TrialForce = forceSensorAnalysis(DAQin, trig, sRate, AVar)

St = trig - AVar.preEveLenQ;
Sp = trig + AVar.posEveLenQ - 1;

coll = DAQin(St:Sp,1);
neck = DAQin(St:Sp,2);

[B,A] = butter(4,40/(sRate/2)); %Low-pass filter under 40

coll  = filter(B,A,abs(coll));
neck  = filter(B,A,abs(neck));

%Eventually some conversion of voltage to force

TrialForce = [coll neck];
end

function [Trialf0Norm] = normf0(Trialf0Raw, f0b)

[r, c] = size(Trialf0Raw);

Trialf0Norm = zeros(r, c);
for i = 1:r
    for j = 1:c
        f = Trialf0Raw(i,j);
        Trialf0Norm(i,j) = (1200*log2(f/f0b));
    end
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

function InflaRespRoute = CalcInflationResponse(auAn, meanTrialf0b, meanRunsf0, limits, plotFolder)
%This calculates the shape of the change in f0 in response to the
%perturbation onset

time   = auAn.time;
curExp = auAn.curExp;
tStep  = auAn.tStep;

disp('Calculating Inflation Response Route')

meanPertMicf0 = meanRunsf0(:,1,2);           %Grabbing the mean mic f0 for all perturbed trials
postOnset     = find(time > 0.5); %Trials are centered at 0.5s before inflation. 
[~, ind]      = min(meanPertMicf0);          %Find the lowest point in whole mean signal

timeFram = postOnset(1):ind;
InflaRespRoute = meanPertMicf0(timeFram);

plotpos = [200 400];
plotdim = [600 600];
InflaRespFig = figure('Color',[1 1 1]);
set(InflaRespFig, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

t = 0:tStep:(tStep*(length(timeFram)-1)); %meanPertMicf0(timeFram,1) - meanPertMicf0(postOnset(1));
plot(t, InflaRespRoute, 'blue', 'LineWidth',3)

xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')
title({'Average Acoustic Response to Inflation'; [curExp '   f0: ' num2str(meanTrialf0b) 'Hz']},...
                     'FontSize', 18,...
                     'FontWeight', 'bold')
axis(limits); box off

set(gca, 'FontSize', 16,...
         'FontWeight','bold');

plTitle = [curExp '_Inflation Response Route.png'];
saveFileName = fullfile(plotFolder, plTitle);
export_fig(saveFileName)
end

function lims = identifyLimits(niAn)

%Full Inidividual Trials: Pressure Sensor
lims.pressure   = [0 4 0 5];

%Aligned Pressure Data
lims.pressureAl = [0 3.5 -0.5 5];

%Full Individual Trials: Force Sensors
lims.force      = [0 4 1 5];

%Full trial f0 analysis
%Full Individual Trials: f0 Audio
if ~isempty(niAn.audioMf0_pPP)
    pertTrials = niAn.audioMf0_pPP;
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

%Section Mean Pertrubed Trials: f0 Audio 
if ~isempty(niAn.audioMf0_meanp)
    [~, Imax] = max(niAn.audioMf0_meanp(:,1)); %Max Pert Onset
    upBoundOn = round(niAn.audioMf0_meanp(Imax,1) + niAn.audioMf0_meanp(Imax,2) + 10);
    [~, Imin] = min(niAn.audioMf0_meanp(:,1)); %Min Pert Onset
    lwBoundOn = round(niAn.audioMf0_meanp(Imin,1) - niAn.audioMf0_meanp(Imin,2) - 10);

    [~, Imax] = max(niAn.audioMf0_meanp(:,3)); %Max Pert Offset
    upBoundOf = round(niAn.audioMf0_meanp(Imax,3) + niAn.audioMf0_meanp(Imax,4) + 10);
    [~, Imin] = min(niAn.audioMf0_meanp(:,3)); %Min Pert Offset
    lwBoundOf = round(niAn.audioMf0_meanp(Imin,3) - niAn.audioMf0_meanp(Imin,4) - 10);

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

function res = packResults(niAn, lims)

res.subject = niAn.subject;
res.run     = niAn.run;
res.curSess = niAn.curSess;
res.AudFB   = niAn.AudFB;

res.numTrials     = niAn.numTrial;
res.numContTrials = niAn.numContTrials;
res.numPertTrials = niAn.numPertTrials;
res.contIdx       = niAn.contIdx;
res.pertIdx       = niAn.pertIdx;
res.pertTrig      = niAn.pertTrig;

res.timeS      = niAn.time_DN;
res.sensorP    = niAn.sensorP_p; %Individual Processed perturbed trials. 
res.lagTimeP   = niAn.lagsPres;
res.lagTimePm  = niAn.meanLagTimeP;
res.riseTimeP  = niAn.riseTimeP;
res.riseTimePm = niAn.riseTimePm;
res.OnOfValP   = niAn.OnOfValP;
res.OnOfValPm  = niAn.OnOfValPm;
res.limitsP    = lims.pressure;

res.timeSAl   = niAn.time_Al;
res.sensorPAl = niAn.sensorP_Al;
res.limitsPAl = lims.pressureAl;

res.timeA     = niAn.time_audio;
res.f0b       = niAn.f0b;

res.numContTrialsPP = niAn.numContTrialsPP;
res.numPertTrialsPP = niAn.numPertTrialsPP;
res.pertTrigPP      = niAn.pertTrigPP;

%Full Individual Trials: Mic/Head f0 Trace 
res.audioMf0TrialPert = niAn.audioMf0_pPP;
res.audioMf0TrialCont = niAn.audioMf0_cPP;
res.audioHf0TrialPert = niAn.audioHf0_pPP;
res.audioHf0TrialCont = niAn.audioHf0_cPP;
res.limitsA           = lims.audio;

%Sections Trials: Mic/Head f0
res.secTime          = niAn.secTime;
res.audioMf0SecPert  = niAn.audioMf0_Secp;
res.audioMf0SecCont  = niAn.audioMf0_Secc;
res.audioHf0SecPert  = niAn.audioHf0_Secp;
res.audioHf0SecCont  = niAn.audioHf0_Secc;

%Mean Sectioned Trials: Mic/Head f0 Trace 
res.audioMf0MeanPert = niAn.audioMf0_meanp; % [MeanSigOn 90%CI MeanSigOff 90%CI]
res.audioMf0MeanCont = niAn.audioMf0_meanc;
res.audioHf0MeanPert = niAn.audioHf0_meanp;
res.audioHf0MeanCont = niAn.audioHf0_meanc;
res.limitsAmean      = lims.audioMean;

%Inflation Response
res.respVar   = niAn.respVar;
res.respVarM  = niAn.respVarMean;
res.respVarSD = niAn.respVarSD;
res.InflaStimVar = niAn.InflaStimVar;
end