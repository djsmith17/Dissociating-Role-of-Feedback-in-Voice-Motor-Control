function OfflineAudapterPlay(varargin)
%This scripts loads a previously recorded audio signal and provides a
%Pitch-shift to it in a similiar fashion that happens during online
%testing.
close all

bf0Vis = 1; bSpVis = 0; bPlay = 0;
if ~isempty(fsic(varargin, '--play'))
    bPlay = 1;
end
 
%Experiment Configurations
expParam.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType       = 'Somatosensory Perturbation_Perceptual';
expParam.subject       = 'Pilot4'; %Subject#, Pilot#, null
expParam.run           = 'Run4';
expParam.curTrial      = [];
expParam.curSubCond    = [];
expParam.defaultGender = 'male';
expParam.masking       = 0;
expParam.bVis          = 0;

dirs = sfDirs(expParam.project, expParam.expType);

dirs.saveFileDir = fullfile(dirs.Data, expParam.subject, expParam.run);
dirs.saveResultsDir = fullfile(dirs.Results, expParam.subject, expParam.run);
dirs.saveFileSuffix = 'offlinePSR';

if exist(dirs.saveResultsDir, 'dir') == 0
    mkdir(dirs.saveResultsDir)
end

%Set up OST and PCF Files
expParam.ostFN = fullfile(dirs.Prelim, 'AFPerturbOST.ost'); check_file(expParam.ostFN);
expParam.pcfFN = fullfile(dirs.Prelim, 'AFPerturbPCF.pcf'); check_file(expParam.pcfFN);

%Should return variables of InflaRespRoute and tStep. 
%Recorded from previous experiments
dirs.InflaRespFile = fullfile(dirs.InflaRespFile, expParam.subject, [expParam.subject '_AveInflaResp.mat']);
try
    load(dirs.InflaRespFile);
catch me
    fprintf('\nSubject Data does not exist at %s \n', dirs.InflaRespFile)
end

d = dir([dirs.saveFileDir, '\*.mat']);
fnames = sort_nat({d.name}); 

%Taking the first trial for ease. File out will be 'data'
load(fullfile(dirs.saveFileDir, fnames{1})); 

Mraw  = data.signalIn; 
fs    = data.params.sRate;
p     = getAudapterDefaultParams('male');



%Level of f0 change based on results from 

p.numTrial    = 1;
p.trialLen    = 2;
p.bPitchShift = 1;
p.dScale      = 2;

trialType = zeros(1,p.numTrial);
trialType(1) = 1;
newInd = randperm(p.numTrial);
trialType = trialType(newInd);

%Resample at 48000Hz
Mraw        = resample(Mraw, data.params.sr * data.params.downFact, fs);
%Split the signal into frames
Mraw_frames = makecell(Mraw, data.params.frameLen * data.params.downFact);

% RMSTHRES = data.rms(svEn(1),1);
% write2pcf(p.ostFN, [RMSTHRES 0.10], 'ost')

for ii = 1:p.numTrial
    expParam.curTrial   = ['Trial' num2str(ii)];
    expParam.curSubCond = [expParam.subject expParam.run expParam.curTrial];
    
    audStimP = setPSRLevels(InflaRespRoute, tStep, expParam.ostFN, expParam.pcfFN, 1);
    
    %Set the OST and PCF functions
%     write2pcf(p.pcfFN, timeWarpFiles(ii,:), 'pcf')
    Audapter('ost', expParam.ostFN, 0);
    Audapter('pcf', expParam.pcfFN, 0);
    
%     Audapter('setParam', 'rmsthr', 5e-3, 0);
    AudapterIO('init', p);
    Audapter('reset');

    for n = 1:length(Mraw_frames)
        Audapter('runFrame', Mraw_frames{n});
    end   
    
    data_offline = AudapterIO('getData');
    Mraw_offline = data_offline.signalIn(1:(end-128)); % Microphone
    Hraw_offline = data_offline.signalOut(129:end);    % Headphones
    fs_offline   = round(data_offline.params.sRate);   % Sampling Rate
    pert  = trialType(ii); 
    
    span = find(data_offline.ost_stat > 0);
     
    span = fs*0.5+1; winL = 0.05; pOve = 0.30;
    [plotf0pts, numPoints, f0_baseline] = sampleParser(Mraw_offline, Hraw_offline, span, fs_offline, winL, pOve);
     
    plotf0pts(:,2) = normf0(plotf0pts(:,2), f0_baseline);
    plotf0pts(:,3) = normf0(plotf0pts(:,3), f0_baseline);
    
    if bf0Vis
        limits = [0 0 0 0];
        drawInterTrialf0(plotf0pts, pert)
    end
    
    if bSpVis 
        OST_MULT = 250;
        visSignals(data_offline, fs, OST_MULT, dirs.saveResultsDir)
    end

    if bPlay; soundsc(data_offline.signalOut, fs); end
    
    fileName = fullfile(dirs.saveResultsDir, expParam.curSubCond, '_', dirs.saveFileSuffix, 'MicIn.wav'); 
    audiowrite(fileName, data_offline.signalIn, fs)
    fileName = fullfile(dirs.saveResultsDir, expParam.curSubCond, '_', dirs.saveFileSuffix, 'HeadOut.wav'); 
    audiowrite(fileName, data_offline.signalOut, fs)
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