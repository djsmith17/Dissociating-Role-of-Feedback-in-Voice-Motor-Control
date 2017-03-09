function OfflineAudapterPlay(varargin)
%This scripts loads a previously recorded audio signal and provides a
%Pitch-shift to it in a similiar fashion that happens during online
%testing.

if isempty(varargin)
else
end

%Experiment Configurations
expParam.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType       = 'Auditory Perturbation_Perceptual';
expParam.subject       = 'Pilot7'; %Subject#, Pilot#, null
expParam.run           = 'Run3';
expParam.numTrial      = 4; %Experimental trials = 40
expParam.curTrial      = [];
expParam.curSubCond    = [];
expParam.perCatch      = 1.00;
expParam.gender        = 'female';
expParam.masking       = 0;
expParam.trialLen      = 4; %Seconds
expParam.bf0Vis        = 0;
expParam.bVis          = 0;
expParam.bPlay         = 0;
expParam.stimType      = 1; %1 for stamped, %2 for sinusoid %3 for linear
expParam.offLineTrial  = 37;

dirs = sfDirs(expParam.project);

dirs.RecFileDir  = fullfile(dirs.RecData, expParam.subject, 'offline');
dirs.RecWaveDir  = fullfile(dirs.RecFileDir, 'wavFiles');

dirs.SavFileDir    = fullfile(dirs.SavData, expParam.subject, expParam.run);
dirs.SavResultsDir = fullfile(dirs.Results, expParam.subject, expParam.run);
dirs.saveFileSuffix = '_offlinePSR';

if exist(dirs.RecFileDir, 'dir') == 0
    mkdir(dirs.RecFileDir)
end
if exist(dirs.RecWaveDir, 'dir') == 0
    mkdir(dirs.RecWaveDir)
end
if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end

%Paradigm Configurations
expParam.sRate              = 48000;  % Hardware sampling rate (before downsampling)
expParam.downFact           = 3;
expParam.sRateAnal          = expParam.sRate/expParam.downFact; %Everything get automatically downsampled! So annoying
expParam.frameLen           = 96;  % Before downsampling
expParam.audioInterfaceName = 'MOTU MicroBook'; %'ASIO4ALL' 'Komplete'

%Set up Audapter
Audapter('deviceName', expParam.audioInterfaceName);
Audapter('setParam', 'downFact', expParam.downFact, 0);
Audapter('setParam', 'sRate', expParam.sRateAnal, 0);
Audapter('setParam', 'frameLen', expParam.frameLen / expParam.downFact, 0);
p = getAudapterDefaultParams(expParam.gender);

%Set up Parameters to control NIDAQ and Perturbatron
s = initNIDAQ;
expParam.sRateQ = s.Rate; %save the sampling rate of the NIDAQ

%Set up OST and PCF Files
expParam.ostFN = fullfile(dirs.Prelim, 'AFPerturbOST.ost'); check_file(expParam.ostFN);
expParam.pcfFN = fullfile(dirs.Prelim, 'AFPerturbPCF.pcf'); check_file(expParam.pcfFN);

%Should return variables of InflaRespRoute and tStep. 
%Recorded from previous experiments
dirs.InflaRespFile = fullfile(dirs.SavData, expParam.subject, [expParam.subject '_AveInflaResp.mat']);
try
    load(dirs.InflaRespFile);
catch me
    fprintf('\nSubject Data does not exist at %s \n', dirs.InflaRespFile)
end

[expParam, p] = setAudFeedType(expParam, dirs, p); %Trials with masking or no... 

expParam.trialType = orderTrials(expParam.numTrial, expParam.perCatch); %numTrials, percentCatch

[expParam.sigs, expParam.trigs] = createPerturbSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType);

%Create a negative voltage signal for the force sensors
negVolSrc = zeros(expParam.sRateQ*expParam.trialLen, 1) - 1;
negVolSrc(1) = 0; negVolSrc(end) = 0;

expParam.resPause = 2.0;

%Taking the first trial for ease. File out will be 'data'
d = dir([dirs.SavFileDir, '\*.mat']);
fnames = sort_nat({d.name}); 
load(fullfile(dirs.SavFileDir, fnames{expParam.offLineTrial})); 

Mraw  = data.signalIn; 
fs    = data.params.sRate;

%Resample at 48000Hz
Mraw        = resample(Mraw, data.params.sr * data.params.downFact, fs);
%Split the signal into frames
Mraw_frames = makecell(Mraw, data.params.frameLen * data.params.downFact);

% RMSTHRES = data.rms(svEn(1),1);

for ii = 1:expParam.numTrial
    expParam.curTrial   = ['Trial' num2str(ii)];
    expParam.curSubCond = [expParam.subject expParam.run expParam.curTrial];
    
    audStimP = setPSRLevels(InflaRespRoute, tStep, expParam.ostFN, expParam.pcfFN, expParam.trialType(ii), expParam.trigs(ii,:,1), expParam.stimType);
    
    %Set the OST and PCF functions
    Audapter('ost', expParam.ostFN, 0);
    Audapter('pcf', expParam.pcfFN, 0);
    
    fprintf('Trial %d\n', ii)
    AudapterIO('init', p);
    Audapter('reset');

    for n = 1:length(Mraw_frames)
        Audapter('runFrame', Mraw_frames{n});
    end
    
    NIDAQsig = [expParam.sigs(:,ii) negVolSrc];
    queueOutputData(s, NIDAQsig);
    
    [dataDAQ, time] = s.startForeground;
    
    data_off = svData(expParam, dirs, p, audStimP, dataDAQ);
    
%     Mraw_off = data_off.signalIn(1:(end-128));  % Microphone
%     Hraw_off = data_off.signalOut(129:end);     % Headphones
%     fs_off   = round(data_off.params.sRate);    % Sampling Rate
%     pert  = expParam.trialType(ii); 
%     
%     span = find(data_off.ost_stat > 0);
%      
%     span = fs*0.5+1; winL = 0.05; pOve = 0.30;
%     [plotf0pts, numPoints, f0_baseline] = sampleParser(Mraw_off, Hraw_off, span, fs_off, winL, pOve);
%      
%     plotf0pts(:,2) = normf0(plotf0pts(:,2), f0_baseline);
%     plotf0pts(:,3) = normf0(plotf0pts(:,3), f0_baseline);
%     
%     if expParam.bf0Vis
%         limits = [0 0 0 0];
%         drawInterTrialf0(plotf0pts, pert)
%     end
%     
%     if expParam.bVis 
%         OST_MULT = 500;
%         visSignals(data_off, fs, OST_MULT, dirs.saveResultsDir)
%     end

    if expParam.bPlay; soundsc(data_off.signalOut, fs); end
    
    pause(expParam.resPause)
end
end

function data = svData(expParam, dirs, p, audStimP, dataDAQ)
%Package all the data into something that is useful for analysis

try
    data = AudapterIO('getData');
    
    data.expParam    = expParam; %Experimental Parameters
    data.dirs        = dirs;     %Directories
    data.p           = p;        %Audapter Parameters
    data.audStimP    = audStimP; %auditory stimulus Parameters
    data.DAQin       = dataDAQ;  %NIDAQ recordings ('Force Sensors')
    save(fullfile(dirs.RecFileDir, [expParam.curSubCond dirs.saveFileSuffix]), 'data')

    audiowrite(fullfile(dirs.RecWaveDir,[expParam.curSubCond dirs.saveFileSuffix '_headOut.wav']), data.signalOut, expParam.sRateAnal)
    audiowrite(fullfile(dirs.RecWaveDir,[expParam.curSubCond dirs.saveFileSuffix '_micIn.wav']), data.signalIn, expParam.sRateAnal)
catch
    disp('Audapter decided not to show up today')
    data = [];
    return
end
end

function [plotf0pts, numPoints, f0_baseline] = sampleParser(mic, head, span, fs, win, oL)
%Finds the value of f0 over windows of the signal 
% St = span - fs*0.5;
% Sp = span + fs*0.7 -1;
% 

St = span(1);

mic = mic(St:end);
head = head(St:end);

numSamp     = length(mic);
AnalysisWin = round(fs*win);
starting    = 1;
noverLap    = AnalysisWin*(1 - oL);

evalSteps = starting:noverLap:(numSamp-AnalysisWin);
numPoints = length(evalSteps);

plotf0pts = [];
for ii = 1:numPoints
    startPt  = evalSteps(ii);
    stopPt   = evalSteps(ii) + AnalysisWin - 1;
    middlePt = round(mean([startPt stopPt]));
    timePt   = (middlePt - 1)/fs;
    
    mic_now  = mic(startPt:stopPt);
    head_now = head(startPt:stopPt);
    
    f0_M = calcf0(mic_now,fs);
    f0_H = calcf0(head_now,fs);
    if f0_M < 50 || f0_M > 300
        f0_M = plotf0pts(ii-1,2);
        disp('Hit')
    end
    
    plotf0pts  = cat(1, plotf0pts, [timePt f0_M f0_H]);
end
f0_baseline = mean(plotf0pts(1:14,2));
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
 
function [normf0pts] = normf0(plotf0pts, Fb)

normf0pts = zeros(size(plotf0pts));
for i = 1:length(plotf0pts)

    F = plotf0pts(i);
    normf0pts(i) = (1200*log2(F/Fb));
end
end

function drawInterTrialf0(plotf0pts, pert)
plotpos = [400 400];
plotdim = [1000 600];
InterTrialNHR = figure('Color', [1 1 1]);
set(InterTrialNHR, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

subplot(2,1,1)
plot(plotf0pts(:,1), plotf0pts(:,2))
xlabel('Time (s)')
ylabel('cents (-)')
title('Pitch_Microphone', 'FontSize', 10, 'FontWeight', 'bold')
axis([0 4 -100 100])
box off   

subplot(2,1,2)
plot(plotf0pts(:,1), plotf0pts(:,3))
xlabel('Time (s)')
ylabel('cents (-)')
title('Pitch_Headphones', 'FontSize', 10, 'FontWeight', 'bold')
axis([0 4 -100 100])
box off

suptitle(['Trial Type: ' num2str(pert)])
                
pause(2)
% close all
end

function visSignals(data, fs, OST_MULT, savedResdir)
plotpos = [200 100];
plotdim = [1000 700];
spectComp = figure('Color', [1 1 1]);
set(spectComp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

frameDur = data.params.frameLen / data.params.sr;
tAxis = 0:frameDur:frameDur * (size(data.rms, 1) - 1);

subplot(2,1,1)
show_spectrogram(data.signalIn, fs, 'noFig');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Microphone In')
box off
plot(tAxis, data.ost_stat * OST_MULT, 'k-');
legend({sprintf('OST status * %d', OST_MULT)});

subplot(2,1,2)
show_spectrogram(data.signalOut, fs, 'noFig');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Headphone Out')
box off

suptitle('OFFLINE ''AFA'' TIME WARP SPECTRUM')

plTitle = 'Offline AFA Time Warp Spectrum';
saveFileName = [savedResdir plTitle '.png'];
export_fig(saveFileName)

end