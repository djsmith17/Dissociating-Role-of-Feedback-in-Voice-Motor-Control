function AFPerturb(varargin)
%Pitch-shift Perturbation experiment. This script measures acoustic output 
%from a participant as they have their auditory feedback perturbed.
%Audapter collects and %manages the recorded acoustic data. This 
%specifically uses a pitch-shift that matches the size of the stimulus seen
%in the somatosensory perturbation experiment.

%This script calls the following (8) functions:
%sfDirs.m
%initNIDAQ.m
%setAudFeedType.m
%orderTrials.m
%createPerturbSignal.m
%setPSRLevels.m
%setPerturbVisualFB.m
%updateVisualFeed.m

%This uses the toolbox from MATLAB-Toolboxes
%speechres

if isempty(varargin)
    targRMS = 55; 
else
    targRMS = varargin{1};
end

%Experiment Configurations
expParam.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType       = 'Auditory Perturbation_Perceptual';
expParam.subject       = 'null'; %Subject#, Pilot#, null
expParam.run           = 'Run1';
expParam.numTrial      = 40; %Experimental trials = 40
expParam.curTrial      = [];
expParam.curSubCond    = [];
expParam.perCatch      = 0.25;
expParam.gender        = 'female';
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
expParam.sRateQ = s.Rate; %save the sampling rate of the NIDAQ

%Set up OST and PCF Files
expParam.ostFN = fullfile(dirs.Prelim, 'AFPerturbOST.ost'); check_file(expParam.ostFN);
expParam.pcfFN = fullfile(dirs.Prelim, 'AFPerturbPCF.pcf'); check_file(expParam.pcfFN);

[expParam, p]      = setAudFeedType(expParam, dirs, p); %Trials with masking or no... ;

expParam.trialType = orderTrials(expParam.numTrial, expParam.perCatch); %numTrials, percentCatch

[expParam.sigs, expParam.trigs] = createPerturbSignal(expParam.trialLen, expParam.numTrial, s.Rate, expParam.sRateAnal, expParam.trialType, expParam.expType);

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

expParam.targRMS   = targRMS; %dB
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
    audStimP = setPSRLevels(InflaRespRoute, tStep, expParam.ostFN, expParam.pcfFN, expParam.trialType(ii), expParam.trigs(ii,:,1));
    
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

    %Grab smooth RMS trace from 'data' structure, compare against baseline
    [color, newPos] = updateVisualFeed(anMsr, data.rms(:,1));
    
    set(rec, 'position', newPos);
    set(rec, 'Color', color); set(rec, 'FaceColor', color);
    set(rec, 'Visible', 'on'); 
    set(fbLines, 'Visible', 'on');  
    
    pause(expParam.resPause)
    set(fbLines, 'Visible', 'off');
    set(rec, 'Visible', 'off');
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

function data = svData(expParam, dirs, s, p, audStimP, dataDAQ)
%Package all the data into something that is useful for analysis

try
    data = AudapterIO('getData');
    
    data.expParam    = expParam; %Experimental Parameters
    data.dirs        = dirs;     %Directories
    data.s           = s;        %NIDAQ Object
    data.p           = p;        %Audapter Parameters
    data.audStimP    = audStimP; %auditory stimulus Parameters
    data.DAQin       = dataDAQ;  %NIDAQ recordings ('Force Sensors')
    save(fullfile(dirs.saveFileDir, expParam.curSubCond), 'data')

    audiowrite(fullfile(dirs.saveWaveDir,[expParam.curSubCond '_headOut.wav']), data.signalOut, expParam.sRateAnal)
    audiowrite(fullfile(dirs.saveWaveDir,[expParam.curSubCond '_micIn.wav']), data.signalIn, expParam.sRateAnal)
catch
    disp('Audapter decided not to show up today')
    data = [];
    return
end
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