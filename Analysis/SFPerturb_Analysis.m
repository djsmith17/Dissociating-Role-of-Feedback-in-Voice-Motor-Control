function SFPerturb_Analysis()
%Analyses the microphone data from the somatosensory perturbation
%experiment. Measures the change in f0 over each trial, and each run for a
%given participant. At the end it approximates a general response to
%inflation to be used in the auditory perturbation experiment

clear all; close all;
%Plot Toggles. This could eventually become an input variable
PltTgl.ForceSensor     = 0;
PltTgl.Trial_time      = 0; %Time-series trial plot
PltTgl.Trial_f0        = 0; %Individual Trial change in NHR
PltTgl.aveTrial_f0     = 0; %Average Trial change in NHR, separated by pert type
PltTgl.aveSessTrial_f0 = 1; 
PltTgl.SPaveSessTrial_f0 = 0;

AVar.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
AVar.expTypes     = {'Somatosensory Perturbation_Perceptual', 'Auditory Perturbation_Perceptual'};
AVar.expInd       = 1; %Either 1 or 2
AVar.curExp       = AVar.expTypes{AVar.expInd};
AVar.participants = {'Pilot7'}; %List of multiple participants
AVar.partiInd     = 1;          %Can select multiple subjs if desired.
AVar.runs         = {'Run1', 'Run2', 'Run3', 'Run4'}; 
AVar.runsInd      = [1 2];
AVar.curRecording = [];

dirs = sfDirs(AVar.project, AVar.curExp);

AVar.anaWinLen   = 0.05; %analysis window length in seconds
AVar.anaWinLenP  = [];    %analysis window length in points
AVar.pOverlap    = 0.30;  %window overlap percentage as decimial

AVar.preEveLen = 0.5; %Amount of time in seconds of observation period before event (onset/offset)
AVar.posEveLen = 1.0; %Amount of time in seconds of observation period after event (onset/offset)
AVar.totEveLen = AVar.preEveLen + AVar.posEveLen; %Total length (seconds) of observation time
AVar.preEveLenP = []; %Amount of points of observation period before event (onset/offset)
AVar.posEveLenP = []; %Amount of points of observation period after event (onset/offset)
AVar.totEveLenP = []; %Total length (points) of observation time
AVar.nOverlap   = []; %Number of points between each analysis window starting indice (Changes with Percent of overlap)
AVar.EvalSteps  = []; %Starting indices for each analysis window
AVar.nEvalSteps = []; %Number of analysis windows;

AVar.svInflaRespRoute = 0;

for i = AVar.partiInd 
    allSessionsf0_St   = [];
    allSessionsf0_Sp   = [];
    allSessionsPert    = [];
    counts             = [0 0];
    for j = AVar.runsInd
        AVar.curRecording   = [AVar.participants{i} ' ' AVar.runs{j}]; %Short hand of experiment details
        dirs.saveFileDir    = fullfile(dirs.Data, AVar.participants{i}, AVar.runs{j}); %Where to find data
        dirs.saveResultsDir = fullfile(dirs.Results, AVar.participants{i}, AVar.runs{j}); %Where to save results
 
        if exist(dirs.saveResultsDir, 'dir') == 0
            mkdir(dirs.saveResultsDir)
        end
        
        %Find total number of files 
        d = dir([dirs.saveFileDir, '\*.mat']);
        AVar.fnames = sort_nat({d.name})';       
        
        limits = [0 1.2 -100 100];
        allplotf0pts_St  = [];
        allplotf0pts_Sp  = [];
        pertRecord       = [];
        countP = 0; countC = 0; %Counting the number of saved perturbed/control trials
       
        for k = 1:length(AVar.fnames)
            %open a given Trial and load 'data.mat' structure
            load([dirs.saveFileDir '\' AVar.fnames{k}]);
            
            %Unpack the 'data.mat' structure
            Mraw  = data.signalIn;  % Microphone
            Hraw  = data.signalOut; % Headphones
            fs    = data.params.sRate;          % Sampling Rate
            pert  = data.expParam.trialType;    % List of trial Order
            span  = data.expParam.spans;        % Pregenerated start and stop points for time-alignment with audio data
            mask  = data.expParam.masking;
            DAQin = data.DAQin;
            audProcDel = data.params.frameLen*4;
            
            ostF  = round(resample(data.ost_stat,32,1));
            ostF  = ostF(129:end);
                           
            AVar.anaWinLenP = round(AVar.anaWinLen*fs);      
            AVar.preEveLenP = round(AVar.preEveLen*fs);  %Amount of points of observation period before event (onset/offset)
            AVar.posEveLenP = round(AVar.posEveLen*fs);  %Amount of points of observation period after event (onset/offset)
            AVar.totEveLenP = AVar.preEveLenP + AVar.posEveLenP; %Total length (points) of observation time
            
            %Determine number of analysis windows and window start indices
            %for frequency analysisover specific period
            AVar.nOverlap   = AVar.anaWinLenP*(1 - AVar.pOverlap); %Number of points between each analysis window starting indice (Changes with Percent of overlap)
            AVar.EvalSteps  = 1:AVar.nOverlap:(AVar.totEveLenP-AVar.anaWinLenP); %Starting indices for each analysis window
            AVar.nEvalSteps = length(AVar.EvalSteps); %Number of analysis windows;


            %saveT decides IF to throw away trial. %base it off of mic data (cleaner)  
            [mic, head, saveT, msg] = preProc(Mraw, Hraw, fs, audProcDel);           
                       
            if saveT == 0 %Don't save the trial :(
                fprintf('Session %d Trial %d not saved. %s\n', j, k, msg)
            elseif saveT == 1 %Save the Trial!
                if pert(k) == 1
                    countP = countP + 1;
                else
                    countC = countC + 1;
                end
                
                %Start of Pert
                [plotf0pts_St, Fb_st] = sampleParser(mic, head, span(k,1), fs, AVar);
                %Stop of Pert
                [plotf0pts_Sp, Fb_sp] = sampleParser(mic, head, span(k,1), fs, AVar); %Short fix in span
                
                prePertInd = plotf0pts_St(:,1) < 0.5;
                f0b = mean(plotf0pts_St(prePertInd,2));
                
                plotf0pts_St(:,2) = normf0(plotf0pts_St(:,2), Fb_st); %Coverted to cents and normalized
                
                plotf0pts_Sp(:,2) = normf0(plotf0pts_Sp(:,2), Fb_st); %Keep Starting baseline frequency
                
                fprintf('Session %d Trial %d saved. %d points \n', j, k, AVar.nEvalSteps)              
                allplotf0pts_St  = cat(3, allplotf0pts_St, plotf0pts_St);
                allplotf0pts_Sp  = cat(3, allplotf0pts_Sp, plotf0pts_St);
                pertRecord       = cat(1, pertRecord, pert(k));
               
                if PltTgl.ForceSensor == 1;
                    sRate = 8000;
                    plot_data_DAQ(sRate, span(k,:), DAQin, AVar.curRecording, dirs.saveResultsDir)
                end
                
                if PltTgl.Trial_time == 1; %Raw time-varying signal
                    drawTrial(Mraw, Hraw, fs, span(k,:))
                end
            
                if PltTgl.Trial_f0 == 1 %Individual Trial change in NHR                   
                    drawIntraTrialf0(plotf0pts_St, plotf0pts_Sp, pert(k), limits, AVar.curRecording, k, dirs.saveResultsDir)
                end
            end          
        end
        curCount = [countC countP];
        
        allSessionsf0_St = cat(3, allSessionsf0_St, allplotf0pts_St);
        allSessionsf0_Sp = cat(3, allSessionsf0_Sp, allplotf0pts_Sp);
        allSessionsPert  = cat(1, allSessionsPert, pertRecord);
        counts = counts + curCount;
               
        %Sort trials of a session by pert type and find averages
        [meanf0pts_St] = sortTrials(allplotf0pts_St, pertRecord);
        [meanf0pts_Sp] = sortTrials(allplotf0pts_Sp, pertRecord);

        %Plots!! See start of script for toggles    
        if PltTgl.aveTrial_f0 == 1      
            drawAVEInterTrialf02(meanf0pts_St, meanf0pts_Sp, limits, curCount, mask, AVar.curRecording, plot_dir)
        end
    end
    
    dirs.saveFileDir = fullfile(dirs.Data, AVar.participants{i}, 'RunAve'); %Where to find data
    AVar.curRecording  = [AVar.participants{i} ' All Runs']; %Short hand of experiment details
    
    if exist(dirs.saveFileDir, 'dir') == 0
        mkdir(dirs.saveFileDir)
    end
    
    %Sort trials of all sessions by pert type and find averages
    [meanSessf0_St] = sortTrials(allSessionsf0_St, allSessionsPert); 
    [meanSessf0_Sp] = sortTrials(allSessionsf0_Sp, allSessionsPert);
    
    %Calculate the response to inflation of the collar. To be used in the
    %Auditory Perturbation Experiment. Only need to use the Average of
    %perturbed Trials
    if AVar.svInflaRespRoute == 1
        [InflaRespRoute, tStep] = CalcInflationResponse(meanSessf0_St{2},1);

        dirs.InflaRespFile = fullfile(dirs.InflaRespFile, AVar.participants{i}, [AVar.participants{i} '_AveInflaResp.mat']);
        save(dirs.InflaRespFile, 'InflaRespRoute', 'tStep')
    end

    %Plots!! See start of script for toggles    
    if PltTgl.aveSessTrial_f0 == 1      
        drawAVEInterTrialf02(meanSessf0_St, meanSessf0_Sp, limits, counts, mask, AVar.curRecording, dirs.saveFileDir)
    end
    
    if PltTgl.SPaveSessTrial_f0 == 1 
        limits = [0 1.2 -80 40];
        drawSPAVEInterTrialf02(meanSessf0_St, meanSessf0_Sp, limits, counts, mask, AVar.curRecording, dirs.saveFileDir)
    end
end
end

function [micP, headP, saveT, msg] = preProc(micR, headR, fs, audProcDel)
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

lenSig = length(x);
t = 0:1/fs:(lenSig-1)/fs;

thresh = 0.5;
[B,A] = butter(4,40/(fs/2));
xenv  = filter(B,A,abs(x));

maxPeak = max(xenv);

I    = find(xenv > thresh*maxPeak);
I0   = I(1);
Iend = length(x); 

[B,A]    = butter(4,(2000)/(fs/2));
filtx    = filtfilt(B,A,x); %Low-pass filtered under 2000Hz
filty    = filtfilt(B,A,y); %Low-pass filtered under 2000Hz
 
micP     = filtx; %Do the analysis on the 
headP    = filty; % use the same indices found for mic (less noise) 

if t(I0) > 1.5
    saveT = 0;  
    msg   = 'We need Gordon!!';
else
    saveT = 1;
    msg   = 'Everything is good'; 
end
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

function [plotf0pts, f0_baseline] = sampleParser(mic, head, span, fs, AVar)
%Finds the value of f0 over windows of the signal 
St = span - AVar.preEveLenP; 
Sp = span + AVar.posEveLenP -1;

try
    mic = mic(St:Sp);
    head = head(St:Sp);
catch
    disp('Sp was too long yo!')
    numSamp = length(mic);
    mic = mic(St:numSamp);
    head = head(St:numSamp);
end   

plotf0pts = [];
for ii = 1:AVar.nEvalSteps
    startPt  = AVar.EvalSteps(ii);
    stopPt   = AVar.EvalSteps(ii) + AVar.anaWinLenP - 1;
    middlePt = round(mean([startPt stopPt]));
    timePt   = (middlePt - 1)/fs;
    
    mic_now  = mic(startPt:stopPt);
    
    f0_M = calcf0(mic_now,fs);
    if f0_M < 50 || f0_M > 300
        f0_M = plotf0pts(ii-1,2);
    end
    
    plotf0pts  = cat(1, plotf0pts, [timePt f0_M]);
end
f0_baseline = mean(plotf0pts(1:14,2));
end

function [normf0pts] = normf0(plotf0pts, Fb)

normf0pts = zeros(size(plotf0pts));
for i = 1:length(plotf0pts)

    F = plotf0pts(i);
    normf0pts(i) = (1200*log2(F/Fb));
end
end

function [meanf0pts] = sortTrials(allplotf0pts, pertRecord)
%This function separates the trials by control or catch trials and
%finds the mean f0 trace and 95% Confidence Interval over multiple trials 
%of a type. 

%Sort Trials by type 
pertVals     = unique(pertRecord(:,1));
num_pertVals = length(pertVals);
time         = allplotf0pts(:,1,1);

meanf0pts  = cell(1, num_pertVals);
for i = 1:num_pertVals
    ind = find(pertRecord(:,1) == pertVals(i));
    TrialsofaType_f0  = allplotf0pts(:,2,ind);
    
    micMean_f0   = mean(squeeze(TrialsofaType_f0(:,1,:)),2);
%     headMean_f0  = mean(squeeze(TrialsofaType_f0(:,2,:)),2);
    micSTD_f0    = std(squeeze(TrialsofaType_f0(:,1,:)),0,2);
%     headSTD_f0   = std(squeeze(TrialsofaType_f0(:,2,:)),0,2); 

    numT = length(ind);
    
    SEM_f0 = micSTD_f0/sqrt(numT); % Standard Error
    CIM_f0 = 1.96*SEM_f0; % 95% confidence Interval
                    
%     ts  = tinv([0.025  0.975],numT-1);      % T-Score
%     CIM = micMean + ts*SEM; 
               
%     ts  = tinv([0.025  0.975],numT-1);      % T-Score
%     CIH = headMean + ts*SEM; 
   meanf0pts{i}    = [time micMean_f0 CIM_f0];
end
end

function [InflaRespRoute, tStep] = CalcInflationResponse(meanSessf0_St, show_plt)
%This calculates the shape of the change in f0 in response to the
%perturbation onset

disp('Calculating Inflation Response Route')

InflaResp = meanSessf0_St(:,1:2); %Grabbing time and mean f0 for perturbed trials
chix = find(InflaResp(:,1) > 0.5); %Trials are centered at 0.5s before inflation. 
[low, ind] = min(InflaResp(:,2)); %Find the lowest point

timeFram = chix(1):ind;
InflaRespRoute = InflaResp(timeFram,:);
tStep = InflaRespRoute(2,1) - InflaRespRoute(1,1);

if show_plt
    plotpos = [200 400];
    plotdim = [600 600];
    InflaRespFig = figure('Color',[1 1 1]);
    set(InflaRespFig, 'Position',[plotpos plotdim],'PaperPositionMode','auto')
    
    t = InflaResp(timeFram,1) - InflaResp(chix(1));
    plot(t,InflaRespRoute)
    
    title('Average Acoustic Response to Inflation')
    xlabel('Time (s)');
    ylabel('f0 (cents)');
    box off
    
end

end

function drawTrial(x, y, fs, span)
h = figure();

figpos = [300 200];
figdim = [1000 500];

set(h,'color',[1 1 1],'position',[figpos figdim])
t = 0:1/fs:((length(x)-1)/fs);
dottedStartx = span(1)/fs*ones(1,50);
dottedEndx   = span(2)/fs*ones(1,50);
dottedy = 5*linspace(-1,1,50);

subplot(1,2,1)
plot(t,x)
hold on
plot(dottedStartx, dottedy,'.k')
hold on
plot(dottedEndx, dottedy,'.k')
xlabel('Time (s)','FontSize', 10, 'FontWeight',  'bold')
ylabel('Amplitude (V)','FontSize', 10, 'FontWeight',  'bold')
title('Microphone', 'FontSize', 10)
axis([0 4 -0.5 0.5])
box off
set(gca, 'FontSize', 10,...
        'FontWeight', 'bold')

subplot(1,2,2)
plot(t,y)
hold on
% plot(dottedStartx, dottedy, '.k')
hold on
% plot(dottedEndx, dottedy,'.k')
xlabel('Time (s)','FontSize', 10, 'FontWeight',  'bold')
ylabel('Amplitude (V)','FontSize', 10, 'FontWeight',  'bold')
title('Headphone', 'FontSize', 10)
axis([0 3 -0.5 0.5])
box off
set(gca, 'FontSize', 10,...
        'FontWeight', 'bold')

% suptitle(['Aspiration Noise Gain of ' num2str(pert)])
end

function drawIntraTrialf0(plotf0pts_St, plotf0pts_Sp, pert, limits, curRecording, k, plotFolder)
plotpos = [200 100];
plotdim = [1300 500];
InterTrialNHR = figure('Color', [1 1 1]);
set(InterTrialNHR, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

dottedStartx = [0.5 0.5];
dottedy = [-300 300];

ha = tight_subplot(1,2,[0.1 0.05],[0.1 0.03],[0.05 0.03]);

axes(ha(1))
plot(plotf0pts_St(:,1), plotf0pts_St(:,2))
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)'); ylabel('f0 (Hz)')

title('Onset of Perturbation', 'FontSize', 10, 'FontWeight', 'bold')
axis(limits); box off

axes(ha(2))
plot(plotf0pts_Sp(:,1), plotf0pts_Sp(:,2))
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)'); ylabel('f0 (Hz)')

title('Offset of Perturbation', 'FontSize', 10, 'FontWeight', 'bold')
axis(limits); box off

suptitle([curRecording ' Trial #' num2str(k)])

if pert == 0
    legend('Unperturbed')
else
    legend('Perturbed')
end

plots = {'IntraTrial_f0'};
for i = 1:length(plots)
    plTitle = [curRecording '_' plots{i} '.png'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end            
end

function drawAVEInterTrialf0(meanf0ptsSt, meanf0ptsSp, limits, counts, mask, curRecording, plotFolder)
plotpos = [200 100];
plotdim = [800 700];
AveInterTrialNHR = figure('Color', [1 1 1]);
set(AveInterTrialNHR, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

dottedStartx = [0.5 0.5];
dottedy = [-300 300];

if mask
    masking = 'With Masking Noise';
    smMask  = 'WMN';
else
    masking = 'Without Masking Noise';
    smMask  = 'WoMN';
end

ha = tight_subplot(2,2,[0.1 0.05],[0.07 0.03],[0.07 0.03]);

axes(ha(1))
errorbar(meanf0ptsSt{1}(:,1), meanf0ptsSt{1}(:,2), meanf0ptsSt{1}(:,3), 'blue')
xlabel('Time (s)')
ylabel('f0 (cents)')
title('Onset of Perturbation: Unperturbed Trials', 'FontSize', 10, 'FontWeight', 'bold')
axis(limits); box off
set(gca,'XTickLabel',{'-0.5' '-0.3' '-0.1' '0.1' '0.3' '0.5' '0.7'})

axes(ha(3))
errorbar(meanf0ptsSt{2}(:,1), meanf0ptsSt{2}(:,2), meanf0ptsSt{2}(:,3), 'black')
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)')
ylabel('f0 (cents)')
title('Onset of Perturbation: Perturbed Trials', 'FontSize', 10, 'FontWeight', 'bold')
axis(limits); box off
set(gca,'XTickLabel',{'-0.5' '-0.3' '-0.1' '0.1' '0.3' '0.5' '0.7'})

axes(ha(2))
errorbar(meanf0ptsSp{1}(:,1), meanf0ptsSp{1}(:,2), meanf0ptsSp{1}(:,3), 'blue')
xlabel('Time (s)')
ylabel('f0 (cents)')
title('Offset of Perturbation: Unperturbed Trials', 'FontSize', 10, 'FontWeight', 'bold')
axis(limits); box off
set(gca,'XTickLabel',{'-0.5' '-0.3' '-0.1' '0.1' '0.3' '0.5' '0.7'})

l1 = legend([num2str(counts(1)) ' Trials']); set(l1,'box', 'off','FontSize', 18);

axes(ha(4))
errorbar(meanf0ptsSp{2}(:,1), meanf0ptsSp{2}(:,2), meanf0ptsSp{2}(:,3), 'black')
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)')
ylabel('f0 (cents)')
title('Offset of Perturbation: Perturbed Trials', 'FontSize', 10, 'FontWeight', 'bold')
axis(limits); box off
set(gca,'XTickLabel',{'-0.5' '-0.3' '-0.1' '0.1' '0.3' '0.5' '0.7'})

suptitle([curRecording ' ' masking])
l2 = legend([num2str(counts(2)) ' Trials']); set(l2,'box', 'off','FontSize', 18);



% 
% leg1 = legend(['Control: ' num2str(counts(2)) ' Trials'], ['Perturbation: ' num2str(counts(1)) ' Trials']);
% set(leg1, 'Position',[0.82,0.88,0.10,0.10]);
                
pause(2)

plots = {'InterTrial_f0'};
for i = 1:length(plots)
    plTitle = [curRecording '_' plots{i} '_' smMask];

    saveFileName = [plotFolder plTitle '.png'];
    export_fig(saveFileName)
end
end

function drawAVEInterTrialf02(meanf0ptsSt, meanf0ptsSp, limits, counts, mask, curRecording, plotFolder)
plotpos = [200 100];
plotdim = [1300 500];
AveInterTrialNHR = figure('Color', [1 1 1]);
set(AveInterTrialNHR, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

dottedStartx = [0.5 0.5];
dottedy = [-300 300];

if mask
    masking = 'With Masking Noise';
    smMask  = 'WMN';
else
    masking = 'Without Masking Noise';
    smMask  = 'WoMN';
end

ha = tight_subplot(1,2,[0.1 0.05],[0.1 0.03],[0.05 0.03]);

axes(ha(1))
errorbar(meanf0ptsSt{1}(:,1), meanf0ptsSt{1}(:,2), meanf0ptsSt{1}(:,3), 'blue') %Unperturbed
hold on
errorbar(meanf0ptsSt{2}(:,1), meanf0ptsSt{2}(:,2), meanf0ptsSt{2}(:,3), 'black') %Perturbed
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)'); ylabel('f0 (cents)')

title('Onset of Perturbation', 'FontSize', 10, 'FontWeight', 'bold')
axis(limits); box off
set(gca,'XTickLabel',{'-0.5' '-0.3' '-0.1' '0.1' '0.3' '0.5' '0.7'})

l0 = legend([num2str(counts(1)) ' Control Trials'], [num2str(counts(2)) ' Perturb Trials']); set(l0,'box', 'off','FontSize', 12);

axes(ha(2))
errorbar(meanf0ptsSp{1}(:,1), meanf0ptsSp{1}(:,2), meanf0ptsSp{1}(:,3), 'blue')  %Unperturbed
hold on
errorbar(meanf0ptsSp{2}(:,1), meanf0ptsSp{2}(:,2), meanf0ptsSp{2}(:,3), 'black') %Perturbed
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)'); ylabel('f0 (cents)')

title('Offset of Perturbation', 'FontSize', 10, 'FontWeight', 'bold')
axis(limits); box off
set(gca,'XTickLabel',{'-0.5' '-0.3' '-0.1' '0.1' '0.3' '0.5' '0.7'})

suptitle([curRecording ' ' masking])
       
pause(2)

plots = {'InterTrial_f0'};
for i = 1:length(plots)
    plTitle = [curRecording '_' plots{i} '_' smMask];

    saveFileName = [plotFolder plTitle '.png'];
    export_fig(saveFileName)
end
end

function drawSPAVEInterTrialf02(meanf0ptsSt, meanf0ptsSp, limits, counts, mask, curRecording, plotFolder)
plotpos = [200 100];
plotdim = [1300 500];
AveInterTrialNHR = figure('Color', [1 1 1]);
set(AveInterTrialNHR, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

dottedStartx = [0.5 0.5];
dottedy = [-300 300];

if mask
    masking = 'With Masking Noise';
    smMask  = 'WMN';
else
    masking = 'Without Masking Noise';
    smMask  = 'WoMN';
end

ha = tight_subplot(1,2,[0.1 0.05],[0.12 0.1],[0.05 0.03]);

axes(ha(1))
errorbar(meanf0ptsSt{1}(:,1), meanf0ptsSt{1}(:,2), meanf0ptsSt{1}(:,3), 'blue', 'LineWidth',2) %Unperturbed
hold on
errorbar(meanf0ptsSt{2}(:,1), meanf0ptsSt{2}(:,2), meanf0ptsSt{2}(:,3), 'black', 'LineWidth',2) %Perturbed
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')

title('Onset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off
set(gca,'XTickLabel',{'-0.5' '-0.3' '-0.1' '0.1' '0.3' '0.5' '0.7'},...
        'FontSize', 16,...
        'FontWeight','bold')

l0 = legend([num2str(counts(1)) ' Control Trials'], [num2str(counts(2)) ' Perturb Trials']); set(l0,'box', 'off','FontSize', 14, 'FontWeight', 'bold');

axes(ha(2))
errorbar(meanf0ptsSp{1}(:,1), meanf0ptsSp{1}(:,2), meanf0ptsSp{1}(:,3), 'blue', 'LineWidth',2)  %Unperturbed
hold on
errorbar(meanf0ptsSp{2}(:,1), meanf0ptsSp{2}(:,2), meanf0ptsSp{2}(:,3), 'black', 'LineWidth',2) %Perturbed
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')

title('Offset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off
set(gca,'XTickLabel',{'-0.5' '-0.3' '-0.1' '0.1' '0.3' '0.5' '0.7'},...
        'FontSize', 16,...
        'FontWeight','bold');

% suptitle([' '])

plots = {'InterTrial_f0_pub'};
for i = 1:length(plots)
    plTitle = [curRecording '_' plots{i} '_' smMask];

    saveFileName = [plotFolder plTitle '.png'];
    export_fig(saveFileName)
end
end

function plot_data_DAQ(sRate, spans, DAQin, curRecording, saveResultsDir)
spans = spans*8000/16000;

[r, c] = size(spans);
pts = length(DAQin);
time = 0:1/sRate:(pts-1)/sRate;

plotpos = [500 500];
plotdim = [1000 400];

sv2File = 0;

for ii = 1:r
    ForceSensorV(ii) = figure('Color', [1 1 1]);
    set(ForceSensorV(ii), 'Position',[plotpos plotdim],'PaperPositionMode','auto')
    
    perturb = zeros(1, pts);
    perturb(spans(ii,1):spans(ii,2)) = -0.5;
    
    subplot(1,2,1)
    plot(time, perturb, 'k')
    hold on
    plot(time, DAQin(:,1,ii), 'b')
    
    xlabel('Time (s)', 'FontSize', 10, 'FontWeight', 'bold') 
    ylabel('Voltage (V)', 'FontSize', 10, 'FontWeight', 'bold')
    title('Collar Sensor', 'FontSize', 10, 'FontWeight', 'bold')
    axis([0 4 -5 5]); box off
    
    subplot(1,2,2)
    plot(time, perturb, 'k')
    hold on
    plot(time, DAQin(:,2,ii), 'b')
    
    xlabel('Time (s)', 'FontSize', 10, 'FontWeight', 'bold')
    ylabel('Voltage (V)', 'FontSize', 10, 'FontWeight', 'bold')
    title('Neck Sensor', 'FontSize', 10, 'FontWeight', 'bold')
    axis([0 4 -5 5]); box off
    
    suptitle('Voltage Change in Force Sensors due to Balloon Inflation')
    
    pltlgd = legend('Perturbation', 'Voltage from Force Sensor');
    set(pltlgd, 'box', 'off',...
                'location', 'best');
   
    set(gca, 'FontSize', 12,...
             'FontWeight', 'bold')
   
    if sv2File == 1
        plTitle = [curRecording  '_ForceSensor_Test ' num2str(ii)];     
        saveFileName = [saveResultsDir plTitle '.png'];
        export_fig(saveFileName) 
    end
    pause()
    close all
end
end