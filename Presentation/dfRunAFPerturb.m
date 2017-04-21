function dfRunAFPerturb(varargin)
%Pitch-shift Perturbation experiment. This script measures acoustic output 
%from a participant as they have their auditory feedback perturbed.
%Audapter collects and %manages the recorded acoustic data. This 
%specifically uses a pitch-shift that matches the size of the stimulus seen
%in the somatosensory perturbation experiment.

%This script calls the following (8) functions:
%dfDirs.m
%initNIDAQ.m
%dfSetAudFB.m
%dfTrialOrder.m
%dfMakePertSignal.m
%dfSetAudapFiles.m
%dfSetVisFB.m
%dfUpdateVisFB.m

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
expParam.run           = 'Run3';
expParam.curExp        = [expParam.subject expParam.run];
expParam.numTrial      = 40; %Experimental trials = 40
expParam.curTrial      = [];
expParam.curSubCond    = [];
expParam.perCatch      = 0.25;
expParam.gender        = 'male';
expParam.masking       = 0;
expParam.trialLen      = 4; %Seconds
expParam.bVis          = 0;
expParam.stimType      = 1; %1 for stamped, %2 for sinusoid %3 for linear

dirs = dfDirs(expParam.project);

dirs.RecFileDir = fullfile(dirs.RecData, expParam.subject, expParam.run);
dirs.RecWaveDir = fullfile(dirs.RecFileDir, 'wavFiles');

if exist(dirs.RecFileDir, 'dir') == 0
    mkdir(dirs.RecFileDir)
end
if exist(dirs.RecWaveDir, 'dir') == 0
    mkdir(dirs.RecWaveDir)
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
[s, niCh, nVS]  = initNIDAQ(expParam.trialLen, 'Dev3');
expParam.sRateQ = s.Rate; % NIDAQ sampling rate
expParam.niCh   = niCh;   % Structure of Channel Names

%Set up OST and PCF Files
expParam.ostFN = fullfile(dirs.Prelim, 'AFPerturbOST.ost'); check_file(expParam.ostFN);
expParam.pcfFN = fullfile(dirs.Prelim, 'AFPerturbPCF.pcf'); check_file(expParam.pcfFN);

[expParam, p]      = dfSetAudFB(expParam, dirs, p); %Trials with masking or no... ;

expParam.trialType = dfSetTrialOrder(expParam.numTrial, expParam.perCatch); %numTrials, percentCatch

[expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType);

%Should give variable of InflaRespRoute. Recorded from previous
%experimentation
dirs.InflaRespFile = fullfile(dirs.SavData, expParam.subject, [expParam.subject '_AveInflaResp.mat']);
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
[anMsr, H1, H2, fbLines, rec, trigCirc] = dfSetVisFB(expParam.targRMS, expParam.boundsRMS, expParam.win);
pause()

DAQin   = [];
rawData = [];
%Close the curtains
pause(1.0) %Let them breathe a sec
for ii = 1:expParam.numTrial
    expParam.curTrial   = ['Trial' num2str(ii)];
    expParam.curExpTrial = [expParam.subject expParam.run expParam.curTrial];
    
    %Level of f0 change based on results from 
    audStimP = dfSetAudapFiles(InflaRespRoute, tStep, expParam.ostFN, expParam.pcfFN, expParam.trialType(ii), expParam.trigs(ii,:,1), expParam.stimType);
    
    %Set the OST and PCF functions
    Audapter('ost', expParam.ostFN, 0);
    Audapter('pcf', expParam.pcfFN, 0);
    
    %Setup which perturb file we want
    NIDAQsig = [expParam.sigs(:,ii) nVS];
    queueOutputData(s, NIDAQsig);
    
    %Cue to begin trial
    set(H1,'Visible','on');
    pause(expParam.cuePause)
    
    %Phonation Start
    set(H1,'Visible','off');
    set(trigCirc,'Visible','on');
    set(H2,'Visible','on');
   
    fprintf('Running Trial %d\n',ii)
    AudapterIO('init', p);
    Audapter('reset');
    Audapter('start');
    
    %This will hold the script for as long as OutputData vector lasts.
    [dataDAQ, time] = s.startForeground;
    
    %Phonation End
    Audapter('stop');
    set(trigCirc,'Visible','off');
    set(H2,'Visible','off'); 
    
    %Save the data
    data = dfSaveRawData(expParam, dirs);
    DAQin = cat(3, DAQin, dataDAQ);
    rawData = cat(1, rawData, data);
    
    %Grab smooth RMS trace from 'data' structure, compare against baseline
    [color, newPos] = dfUpdateVisFB(anMsr, data.rms(:,1));
    
    set(rec, 'position', newPos);
    set(rec, 'Color', color); set(rec, 'FaceColor', color);
    set(rec, 'Visible', 'on'); 
    set(fbLines, 'Visible', 'on');  
    
    pause(expParam.resPause)
    set(fbLines, 'Visible', 'off');
    set(rec, 'Visible', 'off');
end
close all

DRF.dirs        = dirs;
DRF.expParam    = expParam;
DRF.p           = p;
DRF.audStimP    = audStimP;
DRF.DAQin       = DAQin;
DRF.rawData     = rawData; 

dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.curExp dirs.saveFileSuffix '.mat']);
save(dirs.RecFileDir, 'DRF')

if expParam.bVis == 1
    OST_MULT = 500; %Scale factor for OST
    visSignals(data, 16000, OST_MULT, savedWavdir)
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