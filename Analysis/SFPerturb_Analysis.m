function SFPerturb_Analysis()
%Analyses the microphone data from the somatosensory perturbation
%experiment. Measures the change in f0 over each trial, and each run for a
%given participant. At the end it approximates a general response to
%inflation to be used in the auditory perturbation experiment

%12/21/2016: Commented more of the code so I know what it actually does. I
%also made the Inflation Response calculation its own function.
%01/18/2017: Made the function do analysis for the Auditory Perturbation Experiment 
%01/19/2017: Created the variable AVar to localize some important Analysis
%Variables to one structure. This should make passing variables and reading
%values easy. 
%01/24/2017: Changed the name from SFPerturb_Analysis2 to
%LarynPerturb_Analysis. Deleted the file named SFPerturb_Analysis

clear all; close all;
%Plot Toggles. This could eventually become an input variable
PltTgl.Trial_time      = 0; %Time-series trial plot
PltTgl.Trial_f0        = 0; %Individual Trial change in NHR
PltTgl.aveTrial_f0     = 0; %Average Trial change in NHR, separated by pert type
PltTgl.aveSessTrial_f0 = 1; 
PltTgl.SPaveSessTrial_f0 = 1;

AVar.data_folder = 'C:\Users\djsmith\Documents\Pilot Data\Dissociating-Role-of-Feedback-in-Voice-Motor-Control\';
AVar.plot_folder = 'C:\Users\djsmith\Documents\Pilot Results\';

experiments = {'Somatosensory Perturbation_Perceptual', 'Auditory Perturbation_Perceptual'};
exp         = 1; %
participant = {'Pilot6'};
subjs  = 1;             %Can select multiple subjs if desired.
run    = {'Run1', 'Run2', 'Run3', 'Run4'}; 
sess   = [1 2];         %Can select multiple sessions if desired. 

AVar.winLen   = 0.05; %analysis window length in seconds
AVar.pOverlap = 0.30; %Percent Overlap
AVar.anaLen   = 1.20; %What period in time do you want to analyze? Length of vowel is variable.
AVar.nWin     = length(0:AVar.winLen*(1-AVar.pOverlap):(AVar.anaLen-AVar.winLen));

for i = subjs
    allSessionsf0_St   = [];
    allSessionsf0_Sp   = [];
    allSessionsPert    = [];
    counts             = [0 0];
    for j = sess
        AVar.data_dir = [AVar.data_folder '\' experiments{exp} '\' participant{i} '\' run{j} '\']; %Where to find data
        AVar.plot_dir = [AVar.plot_folder '\' experiments{exp} '\' participant{i} '\' run{j} '\']; %Where to save plots
        AVar.Exp_Nm   = [participant{i} ' ' run{j}]; %Short hand of experiment details
        
        if exist(AVar.plot_dir, 'dir') == 0
            mkdir(AVar.plot_dir)
        end
        
        %Find total number of files 
        d = dir([AVar.data_dir, '\*.mat']);
        fnames = sort_nat({d.name})';       
        
        limits = [0 1.2 -100 100];
        allplotf0pts_St  = [];
        allplotf0pts_Sp  = [];
        pertRecord       = [];
        countP = 0; countC = 0; %Counting the number of saved perturbed/control trials
        for k = 1:length(fnames)
            load([AVar.data_dir '\' fnames{k}]);
            Mraw  = data.signalIn(1:(end-128)); % Microphone
            Hraw  = data.signalOut(129:end);    % Headphones
            fs    = round(data.params.sRate);   % Sampling Rate
            pert  = data.trialType;
            span  = data.span/3;
            ostF  = round(resample(data.ost_stat,32,1));
            ostF  = ostF(129:end);
            mask  = data.masking;
            DAQin = data.DAQin;
            
%             St = span(1)-fs*0.5; %0.5s before the start of Pert
%             Sp = span(2)-fs*0.5; %O.5s before the end of Pert
            
            span_m = span - 1.5*fs;
            
%             showMeNIDAQ(DAQin, span, fs)
 
            [mic, head, saveT, msg] = preProc(Mraw, Hraw, fs); %base it off of mic data (cleaner)            
            %saveT decides to throw away the trial or not
                       
            if saveT == 0 %Don't save the trial :(
                fprintf('Session %d Trial %d not saved. %s\n', j, k, msg)
            elseif saveT == 1 %Save the Trial!
                if pert == 1
                    countP = countP + 1;
                else
                    countC = countC + 1;
                end
                
                %Start of Pert
                [plotf0pts_St, numPoints_St, Fb_st] = sampleParser_parts(mic, head, span(1), fs, AVar);
                %Stop of Pert
                [plotf0pts_Sp, numPoints_Sp, Fb_sp] = sampleParser_parts(mic, head, span(2), fs, AVar);
                
                plotf0pts_St(:,2) = normf0(plotf0pts_St(:,2), Fb_st); %Coverted to cents and normalized
                plotf0pts_Sp(:,2) = normf0(plotf0pts_Sp(:,2), Fb_st); %Keep Starting baseline frequency
                
                fprintf('Session %d Trial %d saved. %d points and %d points\n', j, k, numPoints_St, numPoints_Sp)              
                allplotf0pts_St  = cat(3, allplotf0pts_St, plotf0pts_St);
                allplotf0pts_Sp  = cat(3, allplotf0pts_Sp, plotf0pts_Sp);
                pertRecord       = cat(1, pertRecord, pert);
               
                if PltTgl.Trial_time == 1; %Raw time-varying signal
                    drawTrial(Mraw, Hraw, fs, span)
                end
            
                if PltTgl.Trial_f0 == 1 %Individual Trial change in NHR                   
                    drawIntraTrialf0(plotf0pts_St, plotf0pts_Sp, pert, limits, Exp_Nm, k, plot_dir)
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
            drawAVEInterTrialf02(meanf0pts_St, meanf0pts_Sp, limits, curCount, mask, Exp_Nm, plot_dir)
        end
    end
    
    plot_dir = [AVar.plot_folder '\' participant{i} '\RunsAve\'];
    Exp_Nm = [participant{i} ' All Runs'];

    if exist(plot_dir, 'dir') == 0
        mkdir(plot_dir)
    end
    
    %Sort trials of all sessions by pert type and find averages
    [meanSessf0_St] = sortTrials(allSessionsf0_St, allSessionsPert); 
    [meanSessf0_Sp] = sortTrials(allSessionsf0_Sp, allSessionsPert);
    
    %Calculate the response to inflation of the collar. To be used in the
    %Auditory Perturbation Experiment. Only need to use the Average of
    %perturbed Trials
    [InflaRespRoute, tStep] = CalcInflationResponse(meanSessf0_St{2},1);
%     figure; plot(InflaRespRoute)

    InflaRespfile = [AVar.data_folder '\' experiments{exp} '\' participant{subjs} '\' participant{subjs} '_AveInflaResp.mat'];
    save(InflaRespfile, 'InflaRespRoute', 'tStep')

    %Plots!! See start of script for toggles    
    if PltTgl.aveSessTrial_f0 == 1      
        drawAVEInterTrialf02(meanSessf0_St, meanSessf0_Sp, limits, counts, mask, Exp_Nm, plot_dir)
    end
    
    if PltTgl.SPaveSessTrial_f0 == 1 
        limits = [0 1.2 -80 40];
        drawSPAVEInterTrialf02(meanSessf0_St, meanSessf0_Sp, limits, counts, mask, Exp_Nm, plot_dir)
    end
end
end

function showMeNIDAQ(DAQin, span, fs)

t = 0:1/8000:4-1/8000;
span = span/fs;

figure('Color', [1 1 1])
plot(t, DAQin(:,1),'k')
hold on
plot(t,DAQin(:,2),'g')
hold on
plot([span(1) span(1)], [-10 10],'r')
plot([span(2) span(2)], [-10 10],'r')

axis([0 4 -5 5])

box off

end

function [mic1, head1, saveT, msg] = preProc(mic, head, fs)
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
filtx    = filtfilt(B,A,x);
filty    = filtfilt(B,A,y);
 
mic1     = filtx; %Do the analysis on the 
head1    = filty; % use the same indices found for mic (less noise) 

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

function [plotf0pts, numPoints] = sampleParser(mic, head, fs, win, oL)
numSamp    = length(mic);
NHRwindow  = round(fs*win);
starting   = 1;
noverLap   = NHRwindow*(1 - oL);

evalSteps = starting:noverLap:(numSamp-NHRwindow);
numPoints = length(evalSteps);

plotf0pts = [];
for ii = 1:numPoints
    startPt  = evalSteps(ii);
    stopPt   = evalSteps(ii) + NHRwindow - 1;
    middlePt = round(mean([startPt stopPt]));
    timePt   = (middlePt - 1)/fs;
    
    mic_now  = mic(startPt:stopPt);
    head_now = head(startPt:stopPt);
    
    f0_M = calcf0(mic_now,fs);
%     f0_H = calcf0(head_now,fs);
    
    plotf0pts  = cat(1, plotf0pts, [timePt f0_M]);
end
end

function [plotf0pts, numPoints, f0_baseline] = sampleParser_parts(mic, head, span, fs, AVar)
%Finds the value of f0 over windows of the signal 
St = span - fs*0.5; 
Sp = span + fs*0.7 -1;

mic = mic(St:Sp);
head = head(St:Sp);

numSamp    = length(mic);
AnalyisWin = round(fs*AVar.winLen);
starting   = 1;
noverLap   = AnalyisWin*(1 - AVar.pOverlap); %Overlap is in Percent

evalSteps = starting:noverLap:(numSamp-AnalyisWin);
numPoints = length(evalSteps);

plotf0pts = [];
for ii = 1:numPoints
    startPt  = evalSteps(ii);
    stopPt   = evalSteps(ii) + AnalyisWin - 1;
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
axis([0 3 -0.5 0.5])
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
close all
end

function drawIntraTrialf0(plotf0pts_St, plotf0pts_Sp, pert, limits, Exp_Nm, k, plotFolder)
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

suptitle([Exp_Nm ' Trial #' num2str(k)])

if pert == 0
    legend('Unperturbed')
else
    legend('Perturbed')
end


plots = {'IntraTrial_f0'};
for i = 1:length(plots)
    plTitle = [Exp_Nm '_' plots{i}];

    saveFileName = [plotFolder plTitle '.png'];
    export_fig(saveFileName)
end            
end

function drawAVEInterTrialf0(meanf0ptsSt, meanf0ptsSp, limits, counts, mask, Exp_Nm, plotFolder)
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

suptitle([Exp_Nm ' ' masking])
l2 = legend([num2str(counts(2)) ' Trials']); set(l2,'box', 'off','FontSize', 18);



% 
% leg1 = legend(['Control: ' num2str(counts(2)) ' Trials'], ['Perturbation: ' num2str(counts(1)) ' Trials']);
% set(leg1, 'Position',[0.82,0.88,0.10,0.10]);
                
pause(2)

plots = {'InterTrial_f0'};
for i = 1:length(plots)
    plTitle = [Exp_Nm '_' plots{i} '_' smMask];

    saveFileName = [plotFolder plTitle '.png'];
    export_fig(saveFileName)
end
end

function drawAVEInterTrialf02(meanf0ptsSt, meanf0ptsSp, limits, counts, mask, Exp_Nm, plotFolder)
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

suptitle([Exp_Nm ' ' masking])
       
pause(2)

plots = {'InterTrial_f0'};
for i = 1:length(plots)
    plTitle = [Exp_Nm '_' plots{i} '_' smMask];

    saveFileName = [plotFolder plTitle '.png'];
    export_fig(saveFileName)
end
end

function drawSPAVEInterTrialf02(meanf0ptsSt, meanf0ptsSp, limits, counts, mask, Exp_Nm, plotFolder)
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
    plTitle = [Exp_Nm '_' plots{i} '_' smMask];

    saveFileName = [plotFolder plTitle '.png'];
    export_fig(saveFileName)
end
end