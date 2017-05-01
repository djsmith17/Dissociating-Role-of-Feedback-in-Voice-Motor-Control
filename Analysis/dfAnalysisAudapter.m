function [auAn, res] = dfAnalysisAudapter(expParam, rawData, DAQin)
%Analyses the microphone data from the somatosensory perturbation
%experiment. Measures the change in f0 over each trial, and each run for a
%given participant. At the end it approximates a general response to
%inflation to be used in the auditory perturbation experiment

%Require the Signal Processing Toolbox

auAn.curType  = expParam.expType;
auAn.curSubj  = expParam.subject;
auAn.run      = expParam.run;
auAn.curExp   = expParam.curExp; %Short hand of experiment details
auAn.sRate    = expParam.sRateAnal;
auAn.sRateQ   = expParam.sRateQ;
auAn.numTrial = expParam.numTrial;
auAn.trigsT   = expParam.trigs(:,:,1);  %Pregenerated start and stop times for time-alignment with audio data
auAn.trigsA   = expParam.trigs(:,:,3);  %Pregenerated start and stop points (Audapter) for time-alignment with audio data
auAn.trigsQ   = expParam.trigs(:,:,2);  %Pregenerated start and stop points (NIDAQ) for time-alignment with audio data
auAn.trialType = expParam.trialType;    %List of trial Order
auAn.mask      = expParam.masking;

auAn.dnSamp  = 10;
auAn.winLen  = 0.05; %analysis window length in seconds
auAn.pOV     = 0.60;  %window overlap percentage as decimial
auAn.winLenP = auAn.winLen*auAn.sRate;    %analysis window length in points
auAn.tStepP  = auAn.winLenP*(1-auAn.pOV); %Number of points between each analysis window starting indice (Changes with Percent of overlap)
auAn.tStep   = auAn.tStepP/auAn.sRate;            

auAn.actualRecLen = length(rawData(1).signalIn)/auAn.sRate;
auAn.frameT       = linspace(0,auAn.actualRecLen,2053);

auAn.recLen  = expParam.trialLen;
auAn.recLenP = expParam.trialLen*auAn.sRate;

auAn.preEveLen  = 0.5; %Amount of time in seconds of observation period before event (onset/offset)
auAn.posEveLen  = 1.0; %Amount of time in seconds of observation period after event (onset/offset)
auAn.totEveLen  = auAn.preEveLen + auAn.posEveLen; %Total length (seconds) of observation time

auAn.preEveLenP = round(auAn.preEveLen*auAn.sRate);  %Amount of points of observation period before event (onset/offset)
auAn.posEveLenP = round(auAn.posEveLen*auAn.sRate);  %Amount of points of observation period after event (onset/offset)
auAn.totEveLenP = auAn.preEveLenP + auAn.posEveLenP; %Total length (points) of observation time
      
auAn.preEveLenQ = round(auAn.preEveLen*auAn.sRateQ);  %Amount of points of observation period before event (onset/offset) for NIDAQ signal
auAn.posEveLenQ = round(auAn.posEveLen*auAn.sRateQ);  %Amount of points of observation period after event (onset/offset) for NIDAQ signal
auAn.totEveLenQ = auAn.preEveLenQ + auAn.posEveLenQ; %Total length (points_NIDAQ) of observation time
auAn.timeQ      = (0:1:(auAn.totEveLenQ-1))/auAn.sRateQ; %Time points_NIDAQ roughly center of start and stop points of analysis

auAn.winSts  = 1:auAn.tStepP:(auAn.totEveLenP-auAn.winLenP); %Starting indices for each analysis window
auAn.numWin  = length(auAn.winSts); %Number of analysis windows;       

auAn.anaInds(:,1) = auAn.winSts;                      %Start indice for analysis based on EvalStep
auAn.anaInds(:,2) = auAn.winSts + auAn.winLenP - 1;   %Stop indice for analysis based on EvalStep
auAn.time         = mean(auAn.anaInds,2)/auAn.sRate;  %Vector of time points roughly centered on start and stop points of analysis

res.time          = auAn.time;
res.runTrialOrder = [];
res.allTrialf0_St = [];
res.allTrialf0_Sp = [];
res.allTrialf0b   = [];
res.allTrialForce = [];
res.trialCount    = [];
res.meanTrialf0_St = [];
res.meanTrialf0_Sp = [];
res.meanTrialf0b   = [];
res.meanTrialForce_St = [];
res.meanTrialForce_Sp = [];

for ii = 1:auAn.numTrial
    data = rawData(ii);
    
    Mraw = data.signalIn;     % Microphone
    Hraw = data.signalOut;    % Headphones
    OST  = data.ost_stat;
    audProcDel = data.params.frameLen*4;
    
    [mic, head, saveT, saveTmsg] = preProc(Mraw, Hraw, auAn.sRate, audProcDel, auAn.trigsT(ii,1));

    if saveT == 0 %Don't save the trial :(
        fprintf('%s Trial %d not saved. %s\n', auAn.curExp, ii, saveTmsg)
    elseif saveT == 1 %Save the Trial
        fprintf('%s Trial %d saved\n', auAn.curExp, ii)
        
%         Trialf0Raw = signalFrequencyAnalysis(mic, head, auAn.trigsA(ii,1), auAn.sRate, auAn);
        
        
        
        %Start of Pert
        Trialf0Raw_St = signalFrequencyAnalysis(mic, head, auAn.trigsA(ii,1), auAn.sRate, auAn);
        %Stop of Pert
        Trialf0Raw_Sp = signalFrequencyAnalysis(mic, head, auAn.trigsA(ii,2), auAn.sRate, auAn); %When experiment is fixed make this 2!!
        
        trig = auAn.trigsT(ii,1);
        trigSt = trig - 0.5;
        trigSp = trig + 1.0; 
        tFramesSt = find(trigSt >= auAn.frameT); tFrameSt = tFramesSt(end);
        tFramesSp = find(trigSp >= auAn.frameT); tFrameSp = tFramesSp(end);
        frameT_Start = auAn.frameT(tFrameSt:tFrameSp);
        OST_Start = OST(tFrameSt:tFrameSp);
        headedup = find(OST_Start == 3); headUp = headedup(1);
        whenitactuallystarted = frameT_Start(headUp);
        whatTHEDIFF = whenitactuallystarted - trig
        

        prePertInd = auAn.time < 0.5;                    % Grab the first 0.5s, should be no stimulus
        f0b = round(mean(Trialf0Raw_St(prePertInd, 1))); % Baseline fundamental frequency of mic data

        Trialf0Norm_St = normf0(Trialf0Raw_St, f0b); %Coverted to cents and normalized
        Trialf0Norm_Sp = normf0(Trialf0Raw_Sp, f0b); %Coverted to cents and normalized
        
        TrialForce = forceSensorAnalysis(DAQin, auAn.trigsQ(ii,1), auAn.sRateQ, auAn); %At the moment only voltage

        res.runTrialOrder = cat(1, res.runTrialOrder, auAn.trialType(ii));
        res.allTrialf0_St = cat(3, res.allTrialf0_St, Trialf0Norm_St);
        res.allTrialf0_Sp = cat(3, res.allTrialf0_Sp, Trialf0Norm_Sp);
        res.allTrialf0b   = cat(1, res.allTrialf0b, f0b);            %Baseline fundamental frequencies
        res.allTrialForce = cat(3, res.allTrialForce, TrialForce);   %Force sensor values;        
    end
end

%Sort trials within a given run by trial type and find average across trials
[res.meanTrialf0_St, res.meanTrialForce_St, res.trialCount] = sortTrials(res.allTrialf0_St, res.allTrialForce, res.runTrialOrder);
[res.meanTrialf0_Sp, res.meanTrialForce_Sp, res.trialCount] = sortTrials(res.allTrialf0_Sp, res.allTrialForce, res.runTrialOrder);
res.meanTrialf0b = round(mean(res.allTrialf0b,1));

res.f0Limits         = [0 auAn.totEveLen -100 100];
res.InflaRespLimits  = [0 0.5 -200 0];
res.ForceLimits      = [0 auAn.totEveLen 1 3.5];
res.PressureLimits   = [0 auAn.totEveLen 20 30];
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

[B,A]    = butter(4,(2000)/(fs/2));
filtx    = filtfilt(B,A,x); %Low-pass filtered under 2kHz
filty    = filtfilt(B,A,y); %Low-pass filtered under 2kHz
 
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

function f0 = calcf0(x,fs)
% Created by Gabriel Galindo
% Formatted by Dante Smith -12/11/15

lim_inf = ceil(fs/(500));
lim_sup = floor(fs/(50));
U = xcov(x,'unbias');
U = U(ceil(end/2):end);
U = (U(lim_inf:lim_sup)-min(U(lim_inf:lim_sup)))/(max(U(lim_inf:lim_sup)) - min(U(lim_inf:lim_sup)));
[M,P] = findpeaks(U);

if isempty(P)
    f0 = NaN;
else
    P = P(find(M >= 0.9,1,'first'));
    if isempty(P)
        f0 = NaN;
    else
        f0 = fs/(P + lim_inf);
    end

    NFFT = pow2(nextpow2(length(x)/4));
    [Pxx,Fxx] = pwelch(x,NFFT,[],[],fs,'onesided');

    if ~isnan(f0)
        H = Pxx(find(Fxx>=f0,1,'first'));
        if (10*log10(max(Pxx)/H) > 80)
            f0 = NaN;
        end
    end   
end
end

function f0 = calcf02(x, fs)
NFFT = pow2(nextpow2(length(x)/4));

[Pmic, f] = pwelch(x, [], [], [], fs, 'onesided');
    
[val, ind] = max(Pmic);
f0 = f(ind);

end

function Trialf0ResultsRaw = signalFrequencyAnalysis(mic, head, trig, fs, auAn)
%Finds the change in fundamental frequency of windowed signal

%Inputs
%mic:  post-processed single-trial Microphone signal
%head: post-processed signle-trial Headphone signal
%trig: trigger point in audio signals when event occurs (onset/offset)
%fs:   sampling frequency of audio signals
%AVar: structure of analysis variables

%Outputs:
%Trialf0ResultsRaw: Array with rows equal to the number of windows. The first
%column is the time value that the window is centered around. The second
%column is the fundamental frequency of the windowed microphone signal. The
%third column is the fundamental frequency of the windowed headphone
%signal.

St = trig - auAn.preEveLenP; 
Sp = trig + auAn.posEveLenP - 1;

%Grab a big chuck of the signal centered around the event
try
    mic = mic(St:Sp);
    head = head(St:Sp);
catch
    disp('Sp was too long yo!')
    numSamp = length(mic);
    mic = mic(St:numSamp);
    head = head(St:numSamp);
end   

Trialf0ResultsRaw = [];
for ii = 1:auAn.numWin
    startPt  = auAn.anaInds(ii,1);
    stopPt   = auAn.anaInds(ii,2);

    mic_win   = mic(startPt:stopPt);
    head_win  = head(startPt:stopPt);
    
    f0_M = calcf0(mic_win,fs);
    f0_H = calcf0(head_win,fs);
    
%     [f0_time, f0_value, SHR, f0_candidates] = shrp(mic_win, fs);

    
    if f0_M < 50 || f0_M > 300
        fprintf('I calculated a f0 of %d. Replacing it.\n', f0_M)
        f0_M = Trialf0ResultsRaw(ii-1,2);
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

function InflaRespRoute = CalcInflationResponse(AVar, meanTrialf0b, meanRunsf0, limits, plotFolder)
%This calculates the shape of the change in f0 in response to the
%perturbation onset

disp('Calculating Inflation Response Route')

meanPertMicf0 = meanRunsf0(:,1,2);           %Grabbing the mean mic f0 for all perturbed trials
postOnset     = find(AVar.anaTimeVec > 0.5); %Trials are centered at 0.5s before inflation. 
[~, ind]      = min(meanPertMicf0);          %Find the lowest point in whole mean signal

timeFram = postOnset(1):ind;
InflaRespRoute = meanPertMicf0(timeFram);

plotpos = [200 400];
plotdim = [600 600];
InflaRespFig = figure('Color',[1 1 1]);
set(InflaRespFig, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

t = 0:AVar.tStep:(AVar.tStep*(length(timeFram)-1)); %meanPertMicf0(timeFram,1) - meanPertMicf0(postOnset(1));
plot(t, InflaRespRoute, 'blue', 'LineWidth',3)

xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')
title({'Average Acoustic Response to Inflation'; [AVar.curRecording '   f0: ' num2str(meanTrialf0b) 'Hz']},...
                     'FontSize', 18,...
                     'FontWeight', 'bold')
axis(limits); box off

set(gca, 'FontSize', 16,...
         'FontWeight','bold');

plTitle = [AVar.curRecording '_Inflation Response Route.png'];
saveFileName = fullfile(plotFolder, plTitle);
export_fig(saveFileName)
end