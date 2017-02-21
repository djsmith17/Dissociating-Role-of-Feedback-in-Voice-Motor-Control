function AFPerturb(varargin)
%Pitch-shift Perturbation experiment. This specifically uses a pitch-shift
%that matches the size of the stimulus seen in the somatosensory
%perturbation experiment

%12/21/2016: Commented alot of code so I know what the heck I am doing. 
%01/16/2017: Worked towards making the masking noise continuous. I also 
%added more variables to the variables structure for easier access. I also
%made the saving data set of commands its own function. 

%Data Configurations
expParam.subject       = 'null'; %Subject#, Pilot#, null
expParam.run           = 'Run1';
expParam.defaultGender = 'male';
expParam.bVis          = 0;
expParam.masking       = 0;

dirs = sfDirs;

datadir       = 'C:\Users\djsmith\Documents';
expType       = 'Auditory Perturbation_Perceptual'; %%
savedFiledir  = [datadir '\Pilot Data\' expType '\' expParam.subject '\' expParam.run '\'];
savedWavdir   = [savedFiledir '\wavFiles\'];

if exist(savedFiledir, 'dir') == 0
    mkdir(savedFiledir)
end
if exist(savedWavdir, 'dir') == 0
    mkdir(savedWavdir)
end

%Paradigm Configurations
sRate              = 48000;  % Hardware sampling rate (before downsampling)
downFact           = 3;
frameLen           = 96;  % Before downsampling
audioInterfaceName = 'MOTU MicroBook'; %'ASIO4ALL' 'Komplete'

%Set up Audapter
Audapter('deviceName', audioInterfaceName);
Audapter('setParam', 'downFact', downFact, 0);
Audapter('setParam', 'sRate', sRate / downFact, 0);
Audapter('setParam', 'frameLen', frameLen / downFact, 0);
p = getAudapterDefaultParams(expParam.defaultGender);
p.savedFiledir = savedFiledir;
p.savedWavdir  = savedWavdir;
p.postProcSRate= sRate/downFact;

%Set up Parameters to control NIDAQ and Perturbatron
s = initNIDAQ;

%Set up OST and PCF Files
prelimfileLoc = 'C:\Users\djsmith\Documents\MATLAB\Dissociating-Role-of-Feedback-in-Voice-Motor-Control\Presentation\PrelimFiles\';
p.ostFN = [prelimfileLoc 'AFPerturbOST.ost']; check_file(p.ostFN); %%
p.pcfFN = [prelimfileLoc 'AFPerturbPCF.pcf']; check_file(p.pcfFN); %%

%Should give variable of InflaRespRoute. Recorded from previous
%experimentation
InflaRespFile = [datadir '\Pilot Data\Dissociating-Role-of-Feedback-in-Voice-Motor-Control\Somatosensory Perturbation_Perceptual\' expParam.subject '\' expParam.subject '_AveInflaResp.mat'];
try
    load(InflaRespFile);
catch me
    fprintf('\nSubject Data does not exist at %s \n', InflaRespFile)
end

% %Level of f0 change based on results from 
% setPSRLevels(InflaRespRoute, p.ostFN, p.pcfFN);

p.numTrial    = 4; %Experimental trials = 40
p.trialLen    = 4; %Seconds
trialLenPts      = p.trialLen*s.Rate; %seconds converted to points

p = setMasking(p, expParam.masking);

p.trialType = orderTrials(p.numTrial, 0.25); %numTrials, percentCatch

[sigs, spans, spans_t] = createPerturbSignal(s, p.numTrial, trialLenPts, p.trialType, expType);
p.spans = spans*(sRate/s.Rate); %Converting from NIDAQ fs to Audapter fs 

%Create a negative voltage signal for the force sensors
negVolSrc = zeros(s.Rate*p.trialLen, 1) - 1;
negVolSrc(1) = 0; negVolSrc(end) = 0;

%This is where the fun begins
fprintf('\nStarting Trials\n\n')
fprintf('Hit Spacebar when ready\n')
pause()

%Close the curtains
[H1, H2, rec] = createVisualFB();
pause(1.0) %Let them breathe a sec
for ii = 1:p.numTrial
    set(H1,'Visible','on');
    pause(1.0) 
    set(H1,'Visible','off');
    set(H2,'Visible','on');
    
    %Level of f0 change based on results from 
    setPSRLevels(InflaRespRoute, tStep, p.ostFN, p.pcfFN, p.trialType(ii), spans_t(ii,:));
    
    %Set the OST and PCF functions
    Audapter('ost', p.ostFN, 0);
    Audapter('pcf', p.pcfFN, 0);
    
    %Setup which perturb file we want
    NIDAQsig = [sigs(:,ii) negVolSrc];
    queueOutputData(s, NIDAQsig);
   
    AudapterIO('init', p);
    Audapter('reset');
    Audapter('start');
    fprintf('Trial %d\n',ii)
    
    %Play out the Analog Perturbatron Signal. This will hold script for as
    %long as vector lasts. In this case, 4.0 seconds. 
    [data_DAQ, time] = s.startForeground;
    
    Audapter('stop');   
    set(H2,'Visible','off'); 
    
    %Save the data
    data = svData(p, ii, data_DAQ);

    color = chkRMS(data); %How loud were they?    
    set(rec, 'Color', color); set(rec, 'FaceColor', color);
    set(rec, 'Visible','on');  
    
    pause(2.0)
    set(rec, 'Visible','off');
end
close all

if expParam.bVis == 1
    OST_MULT = 500;
    visSignals(data, 16000, OST_MULT, savedWavdir)
end
end

function p = setMasking(p, masking)
%This function sets the type of Auditory Feedback to be played. If it is
%masking noise, then it uses the speech-shaped noise file (SSN.wav) to be
%played.

p.masking = masking;

if masking == 0
    p.fb          = 1;
    p.bPitchShift = 1;
    p.dScale      = 1; %Headphone Scalar
elseif masking == 1
    p.fb          = 2;
%     p.fb3Gain     = 2.0;
    p.bPitchShift = 0;
    p.dScale      = 1; %Headphone Scalar
    noiseWavFN = 'util\SSN.wav'; %Uses Speech-Shaped Noise stored in util
    
    maxPBSize  = Audapter('getMaxPBLen');

    check_file(noiseWavFN);
    [w, fs] = read_audio(noiseWavFN);

    if fs ~= p.sr * p.downFact
        w = resample(w, p.sr * p.downFact, fs);              
    end
    if length(w) > maxPBSize
        w = w(1:maxPBSize);
    end
    Audapter('setParam', 'datapb', w, 1);
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

function color = chkRMS(data)
RMS = mean(data.rms(:,1));

highLim = 100;
lowLim  = 20;

if RMS > highLim
    color = [1 0 0]; %Red, Too loud
elseif RMS < lowLim
    color = [0 0 1]; %Blue, Too soft
else
    color = [0 1 0]; %Green, Just right
end

end

function data = svData(p, ii, data_DAQ)

data = AudapterIO('getData');    
data.ExpVariables  = p;
data.trialType     = p.trialType(ii);
data.span          = p.spans(ii,:);
data.masking       = p.masking;
data.DAQin         = data_DAQ;
save([p.savedFiledir 'Trial' num2str(ii)], 'data')
audiowrite([p.savedWavdir 'Trial' num2str(ii) '_headOut.wav'], data.signalOut, p.postProcSRate)
audiowrite([p.savedWavdir 'Trial' num2str(ii) '_micIn.wav'], data.signalIn, p.postProcSRate)

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