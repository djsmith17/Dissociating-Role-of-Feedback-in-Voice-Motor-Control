function SFPerturb_Analysis()
%Analyses the microphone data from the somatosensory perturbation
%experiment. Measures the change in f0 over each trial, and each run for a
%given participant. At the end it approximates a general response to
%inflation to be used in the auditory perturbation experiment

%Require the Signal Processing Toolbox

clear all; close all; clc
%Plot Toggles. This could eventually become an input variable
PltTgl.ForceSensor     = 0; %Voltage trace of force sensor signal
PltTgl.IntraTrial_T    = 0; %SPL trace of individual trial
PltTgl.IntraTrial_f0   = 0; %f0 trace of individual trial
PltTgl.InterTrial_f0   = 1; %Average f0 trace over all trials of a run
PltTgl.InterRun_f0       = 1; %Average f0 trace over all runs analyzed
PltTgl.InterTrial_AudRes = 1; %Average f0 response trace to auditory pert trials of a run
PltTgl.InterRun_AudRes   = 1; %Average f0 response trace to auditory pert over all runs analyzed
PltTgl.InterTrial_Force  = 0;
PltTgl.InterRun_Force    = 0;
PltTgl.svInflaRespRoute  = 0;

AVar.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
AVar.expTypes     = {'Somatosensory Perturbation_Perceptual', 'Auditory Perturbation_Perceptual'};
AVar.expInd       = 2; %Either 1 or 2
AVar.curExp       = AVar.expTypes{AVar.expInd};
AVar.participants = {'Pilot8'}; %List of multiple participants
AVar.partiInd     = 1;          %Can select multiple subjs if desired.
AVar.runs         = {'Run1', 'Run2', 'Run3', 'Run4', 'offline'}; 
AVar.runsInd      = [3 4];
AVar.curRecording = [];

dirs = sfDirs(AVar.project);
% dirs.saveFileSuffix = '_offlinePSR_Sigmoid';

AVar.anaWinLen   = 0.05; %analysis window length in seconds
AVar.anaWinLenP  = [];    %analysis window length in points
AVar.pOverlap    = 0.60;  %window overlap percentage as decimial

AVar.preEveLen = 0.5; %Amount of time in seconds of observation period before event (onset/offset)
AVar.posEveLen = 1.0; %Amount of time in seconds of observation period after event (onset/offset)
AVar.totEveLen = AVar.preEveLen + AVar.posEveLen; %Total length (seconds) of observation time
AVar.preEveLenP = []; %Amount of points of observation period before event (onset/offset)
AVar.posEveLenP = []; %Amount of points of observation period after event (onset/offset)
AVar.totEveLenP = []; %Total length (points) of observation time
AVar.tStepP     = []; %Number of points between each analysis window starting indice (Changes with Percent of overlap)
AVar.tStep      = []; %time step in seconds between the start of each window 
AVar.EvalSteps  = []; %Starting indices for each analysis window
AVar.nEvalSteps = []; %Number of analysis windows;
AVar.anaInds    = []; %Start and Stop indices for analysis based on EvalStep 
AVar.anaTimeVec = []; %Time points roughly center of start and stop points of analysis
AVar.preEveLenQ = []; %Amount of points of observation period before event (onset/offset) for NIDAQ signal
AVar.posEveLenQ = []; %Amount of points of observation period after event (onset/offset) for NIDAQ signal
AVar.totEveLenQ = []; %Total length (points_NIDAQ) of observation time
AVar.QTimeVec   = []; %Time points_NIDAQ roughly center of start and stop points of analysis

AVar.f0Limits         = [0 AVar.totEveLen -240 180];
AVar.InflaRespLimits  = [0 0.5 -200 0];
AVar.ForceLimits      = [0 AVar.totEveLen 1 3.5];
AVar.PressureLimits   = [0 AVar.totEveLen 20 30];

for i = AVar.partiInd 
    allRunsf0_St   = [];
    allRunsf0_Sp   = [];
    allTrialsOrder = [];
    res.allRunsForce_St = [];
    res.allRunsForce_Sp = [];
    for j = AVar.runsInd
        AVar.curRecording   = [AVar.participants{i} ' ' AVar.runs{j}]; %Short hand of experiment details
        
        dirs.SavFileDir    = fullfile(dirs.SavData, AVar.participants{i}, AVar.runs{j}); %Where to find data
        dirs.SavResultsDir = fullfile(dirs.Results, AVar.participants{i}, AVar.runs{j}); %Where to save results
 
        if exist(dirs.SavResultsDir, 'dir') == 0
            mkdir(dirs.SavResultsDir)
        end
        
        %Find total number of files 
        d = dir([dirs.SavFileDir, '\*.mat']);
        AVar.fnames = sort_nat({d.name})';       
        
        allTrialf0_St  = [];
        allTrialf0_Sp  = [];
        runTrialOrder  = [];
        res.allTrialf0b = [];
        res.allTrialForce = [];
        for k = 1:length(AVar.fnames)
            %open a given Trial and load 'data.mat' structure
            load([dirs.SavFileDir '\' AVar.fnames{k}]);
            
            %Unpack the 'data.mat' structure
            Mraw       = data.signalIn;  % Microphone
            Hraw       = data.signalOut; % Headphones
            fs         = data.params.sRate;          % Sampling Rate
            if isfield(data.expParam, 'sRateQ')
                sRateQ = data.expParam.sRateQ;
            else
                sRateQ = 8000; %from NIDAQ
            end
            
            if isfield(data.expParam, 'spans')
                trigsT      = data.expParam.spansT;  %Pregenerated start and stop times for time-alignment with audio data
                trigsA      = data.expParam.spans;   %Pregenerated start and stop points for time-alignment with audio data
                trigsQ      = trigsA*sRateQ/fs;      %Pregenerated start and stop points (NIDAQ) for time-alignment with audio data
            elseif isfield(data.expParam, 'trigs')
                trigsT      = data.expParam.trigs(:,:,1);  %Pregenerated start and stop times for time-alignment with audio data
                trigsA      = data.expParam.trigs(:,:,3);  %Pregenerated start and stop points (Audapter) for time-alignment with audio data
                trigsQ      = data.expParam.trigs(:,:,2);  %Pregenerated start and stop points (NIDAQ) for time-alignment with audio data
            end
            
            trialType  = data.expParam.trialType;    % List of trial Order
            mask       = data.expParam.masking;
            DAQin      = data.DAQin;
            audProcDel = data.params.frameLen*4;
            
            ostF  = round(resample(data.ost_stat,32,1));
            ostF  = ostF(129:end);
                           
            AVar.anaWinLenP = round(AVar.anaWinLen*fs);      
            AVar.preEveLenP = round(AVar.preEveLen*fs);  %Amount of points of observation period before event (onset/offset)
            AVar.posEveLenP = round(AVar.posEveLen*fs);  %Amount of points of observation period after event (onset/offset)
            AVar.totEveLenP = AVar.preEveLenP + AVar.posEveLenP; %Total length (points) of observation time
            
            %Determine number of analysis windows and window start indices
            %for frequency analysisover specific period
            AVar.tStepP     = AVar.anaWinLenP*(1 - AVar.pOverlap); %Number of points between each analysis window starting indice (Changes with Percent of overlap)
            AVar.tStep      = AVar.tStepP/fs; %time step in seconds between the start of each window
            AVar.EvalSteps  = 1:AVar.tStepP:(AVar.totEveLenP-AVar.anaWinLenP); %Starting indices for each analysis window
            AVar.nEvalSteps = length(AVar.EvalSteps); %Number of analysis windows;
                        
            AVar.anaInds(:,1) = AVar.EvalSteps;                       %Start indice for analysis based on EvalStep 
            AVar.anaInds(:,2) = AVar.EvalSteps + AVar.anaWinLenP - 1; %Stop indice for analysis based on EvalStep
            AVar.anaTimeVec   = mean(AVar.anaInds,2)/fs;              %Vector of time points roughly centered on start and stop points of analysis
            
            AVar.preEveLenQ   = round(AVar.preEveLen*sRateQ);  %Amount of points of observation period before event (onset/offset) for NIDAQ signal
            AVar.posEveLenQ   = round(AVar.posEveLen*sRateQ);  %Amount of points of observation period after event (onset/offset) for NIDAQ signal
            AVar.totEveLenQ   = AVar.preEveLenQ + AVar.posEveLenQ; %Total length (points_NIDAQ) of observation time
            AVar.QTimeVec     = (0:1:(AVar.totEveLenQ-1))/sRateQ; %Time points_NIDAQ roughly center of start and stop points of analysis

%             fprintf('Analysis will be performed over %2.0f bins of length %2.0f points with a %2.0f%% overlap\n', AVar.nEvalSteps, AVar.anaWinLenP, 100*AVar.pOverlap)
            %saveT decides IF to throw away trial. %base it off of mic data (cleaner)  
            [mic, head, saveT, saveTmsg] = preProc(Mraw, Hraw, fs, audProcDel, trigsT(k,1));           
                       
            if saveT == 0 %Don't save the trial :(
                fprintf('Run %d Trial %d not saved. %s\n', j, k, saveTmsg)
            elseif saveT == 1 %Save the Trial!
                
                %Start of Pert
                Trialf0Raw_St = signalFrequencyAnalysis(mic, head, trigsA(k,1), fs, AVar);
                %Stop of Pert
                Trialf0Raw_Sp = signalFrequencyAnalysis(mic, head, trigsA(k,2), fs, AVar); %When experiment is fixed make this 2!!
                
                prePertInd = AVar.anaTimeVec < 0.5;      % Grab the first 0.5s, should be no stimulus
                f0b = round(mean(Trialf0Raw_St(prePertInd, 1))); % Baseline fundamental frequency of mic data
                
                Trialf0Norm_St = normf0(Trialf0Raw_St, f0b); %Coverted to cents and normalized              
                Trialf0Norm_Sp = normf0(Trialf0Raw_Sp, f0b); %Coverted to cents and normalized
                
                fprintf('Run %d Trial %d saved\n', j, k)              
                allTrialf0_St  = cat(3, allTrialf0_St, Trialf0Norm_St);
                allTrialf0_Sp  = cat(3, allTrialf0_Sp, Trialf0Norm_Sp);
                runTrialOrder  = cat(1, runTrialOrder, trialType(k));
                
                res.allTrialf0b   = cat(1, res.allTrialf0b, f0b);    %Baseline fundamental frequencies
                                
                TrialForce = forceSensorAnalysis(DAQin, trigsQ(k,1), sRateQ, AVar); %At the moment only voltage
                res.allTrialForce = cat(3, res.allTrialForce, TrialForce);%Force sensor values;
               
                if PltTgl.ForceSensor == 1; %Voltage trace of force sensor signal
                    drawDAQsignal(sRateQ, trigsQ(k,:), DAQin, AVar.ForceLimits, AVar.curRecording, dirs.SavResultsDir)
                end
                
                if PltTgl.IntraTrial_T == 1; %SPL trace of individual trial
                    drawIntralTrialT(Mraw, Hraw, fs, trigsA(k,:))
                end
            
                if PltTgl.IntraTrial_f0 == 1 %f0 trace of individual trial
                    drawIntraTrialf0(AVar.anaTimeVec, Trialf0Norm_St, Trialf0Norm_Sp, AVar.f0Limits, trialType(k), f0b, AVar.curExp, AVar.curRecording, dirs.SavResultsDir)
                end
            end
        end
        
        %Sort trials within a given run by trial type and find averages
        %across trials
        [meanTrialf0_St, meanTrialForce_St, trialCount] = sortTrials(allTrialf0_St, res.allTrialForce, runTrialOrder);
        [meanTrialf0_Sp, meanTrialForce_Sp, trialCount] = sortTrials(allTrialf0_Sp, res.allTrialForce, runTrialOrder);
        meanTrialf0b = round(mean(res.allTrialf0b,1));
        
        allRunsf0_St   = cat(3, allRunsf0_St, allTrialf0_St);
        allRunsf0_Sp   = cat(3, allRunsf0_Sp, allTrialf0_Sp);
        allTrialsOrder = cat(1, allTrialsOrder, runTrialOrder);
        
        res.allRunsForce_St = cat(3, res.allRunsForce_St, res.allTrialForce);
        res.allRunsForce_Sp = cat(3, res.allRunsForce_Sp, res.allTrialForce);
        
        AVar.ForceLimits      = [0 AVar.totEveLen 1 3.5];
           
        if PltTgl.InterTrial_f0 == 1  %Average f0 trace over all trials of a run 
            drawInterTrialf0(AVar.anaTimeVec, meanTrialf0_St, meanTrialf0_Sp, AVar.f0Limits, trialCount, meanTrialf0b, AVar.curExp, AVar.curRecording, dirs.SavResultsDir)
        end
        
        if PltTgl.InterTrial_AudRes == 1  %Average f0 response trace to auditory pert trials of a run 
            wD = length(trialCount);
            drawInterTrialAudResp(AVar.anaTimeVec, meanTrialf0_St(:,:,wD), meanTrialf0_Sp(:,:,wD), AVar.f0Limits, trialCount(wD), meanTrialf0b, AVar.curExp, AVar.curRecording, dirs.SavResultsDir)
        end
        
        if PltTgl.InterTrial_Force == 1
            drawForceSensorSignal(AVar.QTimeVec, meanTrialForce_St, meanTrialForce_Sp, AVar.ForceLimits, trialCount, AVar.curExp, AVar.curRecording, dirs.SavResultsDir)
        end        
    end
    
    %If I decide to analyze more than 1 run at a time. 
    %Saving myself from over analyzing
    if length(AVar.runsInd) > 1
    
        AVar.curRecording   = [AVar.participants{i} ' All ' AVar.curExp(1:3) ' Runs']; %Short hand of experiment details    
        dirs.SavResultsDir = fullfile(dirs.Results, AVar.participants{i}, 'RunsAve'); %Where to save results
 
        if exist(dirs.SavResultsDir, 'dir') == 0
            mkdir(dirs.SavResultsDir)
        end

        %Sort trials of all sessions by pert type and find averages
        [meanRunsf0_St, meanRunsForce_St, runsCount] = sortTrials(allRunsf0_St, res.allRunsForce_St, allTrialsOrder); 
        [meanRunsf0_Sp, meanRunsForce_Sp, runsCount] = sortTrials(allRunsf0_Sp, res.allRunsForce_Sp, allTrialsOrder);

        %Calculate the response to inflation of the collar. To be used in the
        %Auditory Perturbation Experiment. Only need to use the Average of
        %perturbed Trials

        if PltTgl.svInflaRespRoute == 1
            InflaRespRoute = CalcInflationResponse(AVar, meanTrialf0b, meanRunsf0_St, AVar.InflaRespLimits, dirs.SavResultsDir);
            tStep = AVar.tStep;
            dirs.InflaRespFile = fullfile(dirs.SavData, AVar.participants{i}, [AVar.participants{i} '_AveInflaResp.mat']);
            save(dirs.InflaRespFile, 'InflaRespRoute', 'tStep')
        end

        if PltTgl.InterRun_f0 == 1 %Average f0 trace over all runs analyzed
            drawInterTrialf0(AVar.anaTimeVec, meanRunsf0_St, meanRunsf0_Sp, AVar.f0Limits, runsCount, meanTrialf0b, AVar.curExp, AVar.curRecording, dirs.SavResultsDir)
        end
        
        if PltTgl.InterRun_AudRes == 1 %Average f0 response trace to auditory pert over all runs analyzed
            wD = length(runsCount);
            drawInterTrialAudResp(AVar.anaTimeVec, meanRunsf0_St(:,:,wD), meanRunsf0_Sp(:,:,wD), AVar.f0Limits, runsCount(wD), meanTrialf0b, AVar.curExp, AVar.curRecording, dirs.SavResultsDir)
        end
        
        if PltTgl.InterRun_Force == 1
            drawForceSensorSignal(AVar.QTimeVec, meanRunsForce_St, meanRunsForce_Sp, AVar.ForceLimits, runsCount, AVar.curExp, AVar.curRecording, dirs.SavResultsDir)
        end
    end
end
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

function Trialf0ResultsRaw = signalFrequencyAnalysis(mic, head, trig, fs, AVar)
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

St = trig - AVar.preEveLenP; 
Sp = trig + AVar.posEveLenP - 1;

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
for ii = 1:AVar.nEvalSteps
    startPt  = AVar.anaInds(ii,1);
    stopPt   = AVar.anaInds(ii,2);

    mic_win   = mic(startPt:stopPt);
    head_win  = head(startPt:stopPt);
    
    f0_M = calcf0(mic_win,fs);
    f0_H = calcf0(head_win,fs);
    
    if f0_M < 50 || f0_M > 300
        disp('I had some difficulty calculating f0_M')
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