function dfRunSFPerturb(varargin)
%Laryngeal Perturbation experiment. This script records acoustic output 
%from a participant as they have their larynx physically displaced.
%NIDAQ signal provides Pertrubatron stimulus and Audapter collects and
%manages the recorded acoustic data.

%This script calls the following (7) functions:
%dfDirs.m
%initNIDAQ.m
%dfSetAudFB.m
%dfSetTrialOrder.m
%dfMakePertSignal.m
%dfSetVisFB.m
%dfupdateVisFB.m

%This uses the toolbox from MATLAB-Toolboxes
%speechres

close all;
ET = tic;
rng('shuffle');

if isempty(varargin)
    targRMS = 55; 
else
    targRMS = varargin{1};
end

prompt = {'Subject ID:',...
          'Session ID:',...
          'Gender ("male" or "female")'};
name = 'Subject Information';
numlines = 1;
defaultanswer = {'null', 'SF1', 'female'};
answer = inputdlg(prompt, name, numlines, defaultanswer);

if isempty(answer)
    return
end

num_trials = questdlg('Practice or Full?','Length','Practice','Full','Full') ;
switch num_trials
    case 'Practice'
        numTrials = 4;
        perCatch  = 0.25;
    case 'Full'
        numTrials = 40;
        perCatch  = 0.25;
end

%Experiment Configurations
expParam.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType       = 'Somatosensory Perturbation_Perceptual';
expParam.subject       = answer{1};
expParam.run           = answer{2};
expParam.gender        = answer{3};
expParam.curSess       = [expParam.subject expParam.run];
expParam.numTrial      = numTrials;
expParam.curTrial      = [];
expParam.perCatch      = perCatch;
expParam.masking       = 1;
expParam.trialLen      = 4; %Seconds
expParam.bVis          = 0;
expParam.stimType      = 1; %Always 1. Mirroring the AFPerturb

dirs = dfDirs(expParam.project);
% Folder paths to save data files
dirs.RecFileDir  = fullfile(dirs.RecData, expParam.subject, expParam.run);
dirs.RecWaveDir  = fullfile(dirs.RecFileDir, 'wavFiles');

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
[s, niCh, nVS]  = initNIDAQ(expParam.trialLen, 'Dev2');
expParam.sRateQ = s.Rate; % NIDAQ sampling rate
expParam.niCh   = niCh;   % Structure of Channel Names

%Set up OST and PCF Files
expParam.ostFN = fullfile(dirs.Prelim, 'SFPerturbOST.ost'); check_file(expParam.ostFN);
expParam.pcfFN = fullfile(dirs.Prelim, 'SFPerturbPCF.pcf'); check_file(expParam.pcfFN);

[expParam, p]      = dfSetAudFB(expParam, dirs, p); %Trials with masking or no...  

expParam.trialType = dfSetTrialOrder(expParam.numTrial, expParam.perCatch); %numTrials, percentCatch

[expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, expParam.expType);

expParam.cuePause = 1.0;
expParam.resPause = 2.0;

expParam.targRMS   = targRMS; %dB
expParam.boundsRMS = 3;  %+/- dB
expParam.win       = 2;  %which monitor? 1 or 2

%This is where the fun begins
fprintf('\nStarting Trials\n\n')

%Close the curtains
[anMsr, H1, H2, H3, fbLines, rec, trigCirc] = dfSetVisFB(expParam.targRMS, expParam.boundsRMS, expParam.win);

%Close the curtains
pause(5); %Let them breathe a sec
set(H3,'Visible','off');

DAQin = []; rawData = [];
for ii = 1:expParam.numTrial
    expParam.curTrial     = ['Trial' num2str(ii)];
    expParam.curSessTrial = [expParam.subject expParam.run expParam.curTrial];
    
    %Used later in audio version
    audStimP = [];
        
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
    
    fprintf('Trial %d\n',ii)
    AudapterIO('init', p);
    Audapter('reset');
    Audapter('start');
    
    %Play out the Analog Perturbatron Signal. This will hold script for as
    %long as vector lasts. In this case, 4.0 seconds. 
    [dataDAQ, time] = s.startForeground;
     
    %Phonation End
    Audapter('stop');
    set(trigCirc,'Visible','off');
    set(H2,'Visible','off');
    
    %Save the data
    data = dfSaveRawData(expParam, dirs);
    DAQin   = cat(3, DAQin, dataDAQ);
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
close all;
elapsed_time = toc(ET)/60;
fprintf('\nElapsed Time: %f (min)\n', elapsed_time)

expParam.elapsedTime = elapsed_time;
DRF.dirs        = dirs;
DRF.expParam    = expParam;
DRF.p           = p;
DRF.audStimP    = audStimP;
DRF.DAQin       = DAQin;
DRF.rawData     = rawData; 

dirs.RecFileDir = fullfile(dirs.RecFileDir, [expParam.subject expParam.run dirs.saveFileSuffix 'DRF.mat']);
switch num_trials
    case 'Practice'
        return
    case 'Full'
        save(dirs.RecFileDir, 'DRF'); %Only save if it was a full set of trials
end

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