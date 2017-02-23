function AFPerturb(varargin)
%Pitch-shift Perturbation experiment. This script measures acoustic output 
%from a participant as they have their auditory feedback perturbed.
%Audapter collects and %manages the recorded acoustic data. This 
%specifically uses a pitch-shift that matches the size of the stimulus seen
%in the somatosensory perturbation experiment.

%This calls the functions:
%sfDirs.m
%initNIDAQ.m
%setAudFeedType.m
%orderTrials.m
%createPerturbSignal.m
%setPSRLevels.m

%This uses the toolbox from MATLAB-Toolboxes
%speechres

%Experiment Configurations
expParam.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType       = 'Auditory Perturbation_Perceptual';
expParam.subject       = 'null'; %Subject#, Pilot#, null
expParam.run           = 'Run1';
expParam.numTrial      = 40; %Experimental trials = 40
expParam.curTrial      = [];
expParam.curSubCond    = [];
expParam.perCatch      = 0.25;
expParam.gender        = 'male';
expParam.masking       = 0;
expParam.trialLen      = 4; %Seconds
expParam.bVis          = 0;

dirs = sfDirs(expParam.project, expParam.expType);

dirs.saveFileDir = fullfile(dirs.Data, expParam.subject, expParam.run);
dirs.saveWaveDir = fullfile(dirs.saveFileDir, 'wavFiles');

if exist(dirs.saveFileDir, 'dir') == 0
    mkdir(dirs.saveFileDir)
end
if exist(dirs.saveWaveDir, 'dir') == 0
    mkdir(dirs.saveWaveDir)
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

%Set up OST and PCF Files
expParam.ostFN = fullfile(dirs.Prelim, 'AFPerturbOST.ost'); check_file(expParam.ostFN);
expParam.pcfFN = fullfile(dirs.Prelim, 'AFPerturbPCF.pcf'); check_file(expParam.pcfFN);

[expParam, p]      = setAudFeedType(expParam, dirs, p); %Trials with masking or no... ;

expParam.trialType = orderTrials(expParam.numTrial, expParam.perCatch); %numTrials, percentCatch

[expParam.sigs, expParam.spans, expParam.spansT] = createPerturbSignal(expParam.trialLen, expParam.numTrial, s.Rate, expParam.trialType, expParam.expType);
expParam.spans = expParam.spans*(expParam.sRateAnal / s.Rate); %Converting from NIDAQ fs to Audapter analysis fs 

%Create a negative voltage signal for the force sensors
negVolSrc = zeros(s.Rate*expParam.trialLen, 1) - 1;
negVolSrc(1) = 0; negVolSrc(end) = 0;

%Should give variable of InflaRespRoute. Recorded from previous
%experimentation
dirs.InflaRespFile = fullfile(dirs.InflaRespFile, expParam.subject, [expParam.subject '_AveInflaResp.mat']);
try
    load(dirs.InflaRespFile);
catch me
    fprintf('\nSubject Data does not exist at %s \n', dirs.InflaRespFile)
end

expParam.cuePause = 1.0;
expParam.resPause = 2.0;

expParam.targRMS   = 55; %***dB Example at the moment***
expParam.boundsRMS = 3;  %+/- dB
expParam.win       = 2;  %which monitor? 1 or 2

%This is where the fun begins
fprintf('\nStarting Trials\n\n')
fprintf('Hit Spacebar when ready\n')

%Close the curtains
[anMsr, H1, H2, fbLines, rec] = setPerturbVisualFB(expParam.targRMS, expParam.boundsRMS, expParam.win);
pause()

%Close the curtains
pause(1.0) %Let them breathe a sec
for ii = 1:expParam.numTrial
    expParam.curTrial   = ['Trial' num2str(ii)];
    expParam.curSubCond = [expParam.subject expParam.run expParam.curTrial];
    
    %Level of f0 change based on results from 
    audStimP = setPSRLevels(InflaRespRoute, tStep, expParam.ostFN, expParam.pcfFN, expParam.trialType(ii), expParam.spansT(ii,:));
    
    %Set the OST and PCF functions
    Audapter('ost', expParam.ostFN, 0);
    Audapter('pcf', expParam.pcfFN, 0);
    
    %Setup which perturb file we want
    NIDAQsig = [expParam.sigs(:,ii) negVolSrc];
    queueOutputData(s, NIDAQsig);
    
    %Cue to begin trial
    set(H1,'Visible','on');
    pause(expParam.cuePause)
    
    %Phonation
    set(H1,'Visible','off');
    set(H2,'Visible','on');
   
    fprintf('Trial %d\n',ii)
    AudapterIO('init', p);
    Audapter('reset');
    Audapter('start');
    
    %Play out the Analog Perturbatron Signal. This will hold script for as
    %long as vector lasts. In this case, 4.0 seconds. 
    [dataDAQ, time] = s.startForeground;
    
    Audapter('stop');   
    set(H2,'Visible','off'); 
    
    %Save the data
    data = svData(expParam, dirs, s, p, audStimP, dataDAQ);

    [color, newPos] = updateVisualFeed(anMsr, data.rms);
    
    set(rec, 'position', newPos);
    set(rec, 'Color', color); set(rec, 'FaceColor', color);
    set(rec, 'Visible','on'); 
    set(fbLines, 'Visible', 'on');  
    
    pause(expParam.resPause)
    set(fbLines, 'Visible', 'off');
    set(rec, 'Visible','off');
end
close all

if expParam.bVis == 1
    OST_MULT = 500; %Scale factor for OST
    visSignals(data, 16000, OST_MULT, savedWavdir)
end
end

function plotPerturb(s, lenT, sig)

t = (0:1:lenT-1)/s.Rate;
plot(t, sig);

xlabel('Time'); 
ylabel('Voltage'); 
legend('Analog Output 0');
end

function [H1, H2, rec] = createVisualFB()
% figure0 = figure('NumberTitle','off','Color',[0 0 0],'Position',[0 0 1920 1080],'MenuBar','none');

figure1 = figure('NumberTitle','off','Color',[0 0 0],'Position',[1920 0 1681 1050],'MenuBar','none');

H1 = annotation(figure1,'textbox',[0.46 0.46 0.2 0.2],...
                        'Color',[1 1 1],...
                        'String',{'+'},...
                        'LineStyle','none',...
                        'FontWeight','bold',...
                        'FontSize',160,...
                        'FontName','Arial',...
                        'BackgroundColor',[0 0 0],...
                        'Visible','on');

H2 = annotation(figure1,'textbox',[0.38 0.46 0.2 0.2],...
                        'Color',[1 1 1],...
                        'String',{'eee'},...
                        'LineStyle','none',...
                        'FontWeight','bold',...
                        'FontSize',160,...
                        'FontName','Arial',...
                        'BackgroundColor',[0 0 0],...
                        'visible','off');
                    
rec = annotation(figure1,'rectangle',[0.25 0.1 0.5 0.1],...
                         'Color',[0 1 0],...
                         'LineStyle','none',...
                         'FaceColor',[0 1 0],...
                         'visible','off');

end

function color = chkRMS(anMsr, rms)
RMS = mean(rms(:,1));

%based on data I found the new height to be this
newRecHeight = 0.55;
newDrawHeight = newRecHeight + anMsr.bMar;

newPos = [anMsr.recXSt anMsr.recYSt anMsr.recWidth newRecHeight];

if newDrawHeight > anMsr.drawMaxH
    color = 'red';
elseif newDrawHeight < anMsr.drawMinH
    color = 'red';
else
    color = 'green';
end
end

function data = svData(expParam, dirs, s, p, audStimP, dataDAQ)
%Package all the data into something that is useful for analysis

data = AudapterIO('getData');   
data.expParam    = expParam;
data.dirs        = dirs;
data.s           = s;
data.p           = p;
data.audStimP    = audStimP;
data.DAQin       = dataDAQ;
save(fullfile(dirs.saveFileDir, expParam.curSubCond), 'data')

audiowrite(fullfile(dirs.saveWaveDir,[expParam.curSubCond '_headOut.wav']), data.signalOut, p.postProcSRate)
audiowrite(fullfile(dirs.saveWaveDir,[expParam.curSubCond '_micIn.wav']), data.signalIn, p.postProcSRate)
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

suptitle('ONLINE ''AFA'' TIME WARP SPECTRUM')

plTitle = 'Online AFA Time Warp Spectrum';
saveFileName = [savedResdir plTitle '.png'];

% export_fig(saveFileName)
end