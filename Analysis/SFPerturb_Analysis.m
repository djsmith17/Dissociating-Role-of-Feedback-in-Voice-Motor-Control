function SFPerturb_Analysis()

clear all; close all;
pltTgl.Trial_time      = 0; %Time-series trial plot
pltTgl.Trial_f0        = 0; %Individual Trial change in NHR
pltTgl.aveTrial_f0     = 1; %Average Trial change in NHR, separated by pert type
pltTgl.aveSessTrial_f0 = 0; 

data_folder = 'C:\Users\djsmith\Documents\Pilot Data\Somatosensory Perturbation_Perceptual';
plot_folder = 'C:\Users\djsmith\Documents\Pilot Results\Somatosensory Perturbation_Perceptual';
participant = {'Pilot1','Pilot2'};
subjs  = 1; %Can select multiple subjs if desired.
run    = {'Session1', 'Session2', 'Session3', 'Session4', 'Session5'}; 
sess   = 1; %Can select multiple sessions if desired. 

winL   = 0.05; %analysis window length in seconds
pOve   = 0.30; %Percent Overlap
anaLen = 2.00; %What period in time do you want to analyze? Length of vowel is variable.
nWin   = length(0:winL*(1-pOve):(anaLen-winL));

for i = subjs
    allSessionsf0     = [];
    for j = sess
        data_dir = [data_folder '\' participant{i} '\' run{j} '\']; %Where to find data
        plot_dir = [plot_folder '\' participant{i} '\' run{j} '\']; %Where to save plots
        Exp_Nm   = [participant{i} ' ' run{j}]; %Short hand of experiment details
        
        if exist(plot_dir, 'dir') == 0
            mkdir(plot_dir)
        end
        
        %Find total number of files
        d = dir([data_dir, '\*.mat']);
        fnames = sort_nat({d.name});       
        
        limits = [0 2 -60 60];
        allplotf0pts  = [];
        pertRecord    = [];
        countP = 0; countC = 0;
        for k = 1:length(fnames)
            load([data_dir '\' fnames{k}]);
            Mraw  = data.signalIn(1:(end-128)); % Microphone
            Hraw  = data.signalOut(129:end);    % Headphones
            fs    = round(data.params.sRate);   % Sampling Rate
            pert  = data.trialType;
            span  = data.span/3;
            ostF  = round(resample(data.ost_stat,32,1));
            ostF  = ostF(129:end);
            mask  = data.masking;
            DAQin = data.DAQin;
            
            St = span(1)-fs*0.5; %0.5s before the start of Pert
            Sp = span(2)+fs*0.5; %O.5s after the end of Pert
            
            span_m = span - 0.5*fs;
            
            showMeNIDAQ(DAQin, span, fs)
 
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
                                
                [plotf0pts, numPoints] = sampleParser(mic(St:Sp), head(St:Sp), fs, winL, pOve);
                
                plotf0pts(:,2) = normf0(plotf0pts(:,2)); %Coverted to cents and normalized
                
                fprintf('Session %d Trial %d saved. %d points\n', j, k, numPoints)              
                allplotf0pts  = cat(3, allplotf0pts, plotf0pts);
                pertRecord    = cat(1, pertRecord, pert);
               
                if pltTgl.Trial_time == 1; %Raw time-varying signal
                    drawTrial(Mraw, Hraw, fs, span)
                end
            
                if pltTgl.Trial_f0 == 1 %Individual Trial change in NHR                   
                    drawIntraTrialf0(plotf0pts, pert, fs, span_m, limits, Exp_Nm, plot_dir)
                end
            end          
        end
        allSessionsf0 = cat(3, allSessionsf0, allplotf0pts);
        counts = [countP countC];
               
        %Sort trials of a session by pert type and find averages
        [meanf0pts] = sortNHRTrials(allplotf0pts, pertRecord); 

        %Plots!! See start of script for toggles    
        if pltTgl.aveTrial_f0 == 1      
            drawAVEInterTrialf0(meanf0pts, limits, counts, mask, Exp_Nm, plot_dir)
        end
    end
    
    plot_dir = [plot_folder '\' participant{i} '\SessionsAve\'];
    Exp_Nm = [participant{i} '_AllSessions'];

    if exist(plot_dir, 'dir') == 0
        mkdir(plot_dir)
    end
    
    %Sort trials of all sessions by pert type and find averages
    [meanSessf0] = sortNHRTrials(allSessionsf0, pertRecord); 

    %Plots!! See start of script for toggles    
    if pltTgl.aveSessTrial_f0 == 1      
        drawAVEInterTrialf0(meanSessf0, Exp_Nm, plot_dir)
    end
end
end

function showMeNIDAQ(DAQin, span, fs)

t = 0:1/8000:4-1/8000;
span = span/fs;

figure
plot(t, DAQin(:,1),'k')
hold on
plot(t,DAQin(:,2),'g')
hold on
plot([span(1) span(1)], [-10 10],'r')
plot([span(2) span(2)], [-10 10],'r')

axis([0 4 -5 5])

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

if t(I0) > 0.5
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

function [normf0pts] = normf0(plotf0pts)

Fb = mean(plotf0pts(1:14));

normf0pts = zeros(size(plotf0pts));
for i = 1:length(plotf0pts)

    F = plotf0pts(i);
    normf0pts(i) = (1200*log2(F/Fb));
end
end

function [meanf0pts] = sortNHRTrials(allplotf0pts, pertRecord)
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
    CIM_f0 = 1.96*SEM_f0; 
                    
%     ts  = tinv([0.025  0.975],numT-1);      % T-Score
%     CIM = micMean + ts*SEM; 
               
%     ts  = tinv([0.025  0.975],numT-1);      % T-Score
%     CIH = headMean + ts*SEM; 
   meanf0pts{i}    = [time micMean_f0 CIM_f0];
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

function drawIntraTrialf0(plotf0pts, pert, fs, span, limits, Exp_Nm, plotFolder)
plotpos = [400 400];
plotdim = [1000 500];
InterTrialNHR = figure('Color', [1 1 1]);
set(InterTrialNHR, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

dottedStartx = [span(1) span(1)]/fs;
dottedEndx   = [span(2) span(2)]/fs;
dottedy = [-300 300];

plot(plotf0pts(:,1), plotf0pts(:,2))
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
hold on
plot(dottedEndx, dottedy,'k','LineWidth',4)
xlabel('Time (s)')
ylabel('f0 (Hz)')
title('Pitch', 'FontSize', 10, 'FontWeight', 'bold')
axis(limits)
box off
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

function drawAVEInterTrialf0(meanf0pts, limits, counts, mask, Exp_Nm, plotFolder)
plotpos = [400 200];
plotdim = [900 700];
AveInterTrialNHR = figure('Color', [1 1 1]);
set(AveInterTrialNHR, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

dottedStartx = [0.5 0.5];
dottedEndx   = [1.5 1.5];
dottedy = [-300 300];

if mask
    masking = 'With Masking Noise';
else
    masking = 'Without Masking Noise';
end

subplot(2,1,1)
errorbar(meanf0pts{1}(:,1), meanf0pts{1}(:,2), meanf0pts{1}(:,3), 'blue')
xlabel('Time (s)')
ylabel('f0 (cents)')
title('Pitch: Unperturbed', 'FontSize', 10, 'FontWeight', 'bold')
axis(limits)
box off
suptitle([Exp_Nm ' ' masking])
legend([num2str(counts(2)) ' Trials'])

subplot(2,1,2)
errorbar(meanf0pts{2}(:,1), meanf0pts{2}(:,2), meanf0pts{2}(:,3), 'black')
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
hold on
plot(dottedEndx, dottedy,'k','LineWidth',4)
xlabel('Time (s)')
ylabel('f0 (cents)')
title('Pitch: Perturbed', 'FontSize', 10, 'FontWeight', 'bold')
axis(limits)
box off
legend([num2str(counts(1)) ' Trials'])
% 
% leg1 = legend(['Control: ' num2str(counts(2)) ' Trials'], ['Perturbation: ' num2str(counts(1)) ' Trials']);
% set(leg1, 'Position',[0.82,0.88,0.10,0.10]);
                
pause(2)

plots = {'InterTrial_f0'};
for i = 1:length(plots)
    plTitle = [Exp_Nm '_' plots{i}];

    saveFileName = [plotFolder plTitle '.png'];
    export_fig(saveFileName)
end
end