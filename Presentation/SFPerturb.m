function SFPerturb(varargin)
%This script is used to present the behavioral portion for the laryngeal
%perturbation experiment. This uses Audapter. 

%09/09/2016: Made the masking setup its own function for readability
%12/21/2016: Commented more of the code so I actually know what it does
%01/16/2017: Worked towards making the masking noise continuous. I also 
%added more variables to the variables structure for easier access. I also
%made the saving data set of commands its own function. 
%01/24/2017: Added functionality to take in force sensor data. Also renamed
%the function from SFPerturb3 to SFPerturb.

%Data Configurations
subject       = 'null'; %Subject#, Pilot#, null
run           = 'Run1';
defaultGender = 'male';
bVis          = 0;
masking       = 1;

dirs = sfDirs;

datadir       = 'C:\Users\djsmith\Documents';
expType       = 'Pilot Data\Somatosensory Perturbation_Perceptual';
savedFiledir  = [datadir '\' expType '\' subject '\' run '\'];
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
p = getAudapterDefaultParams(defaultGender);
p.savedFiledir = savedFiledir;
p.savedWavdir  = savedWavdir;
p.postProcSRate= sRate/downFact;

%Set up Perturbatron
s = initPerturb;

%Set up OST and PCF Files
prelimfileLoc = 'C:\Users\djsmith\Documents\MATLAB\SFPerturb\PrelimFiles\';
p.ostFN = [prelimfileLoc 'SFPerturbOST.ost']; check_file(p.ostFN);
p.pcfFN = [prelimfileLoc 'SFPerturbPCF.pcf']; check_file(p.pcfFN);

p.numTrial    = 4; %Experimental trials = 40
p.trialLen    = 4; %Seconds
trialLenPts   = p.trialLen*s.Rate; %seconds converted to points

p = setMasking(p, masking); %Trials with masking or no...  

p.trialType = orderTrialsLaced(p.numTrial, 0.25); %numTrials, percentCatch

[sigs, spans] = createPerturbs(s, p.numTrial, trialLenPts, p.trialType);
p.spans = spans*(sRate/s.Rate); %Converting from NIDAQ fs to Audapter fs 

%This is where the fun begins
fprintf('\nStarting Trials\n\n')
fprintf('Hit Spacebar when ready\n')
[H1, H2, rec] = createVisualFB();
pause()

%Close the curtains
pause(1.0) %Let them breathe a sec
for ii = 1:p.numTrial
    %Set the OST and PCF functions
    Audapter('ost', p.ostFN, 0);
    Audapter('pcf', p.pcfFN, 0);
    
    %Setup which perturb file we want
    queueOutputData(s, sigs(:,ii));
    
    %Cue to begin trial
    set(H1,'Visible','on');
    pause(1.0)
    
    %Phonation
    set(H1,'Visible','off');
    set(H2,'Visible','on');  
    
    AudapterIO('init', p);
    Audapter('reset');
    Audapter('start');
    fprintf('Trial %d\n',ii)
    
    %Play out the Analog Perturbatron Signal. This will hold script for as
    %long as vector lasts. In this case, 4.0 seconds. 
    [data_DAQ, time] = s.startForeground;
     
    Audapter('stop');  
    set(H2,'Visible','off');
    
    data = svData(p, ii, data_DAQ);
    
%     data = AudapterIO('getData');   
%     data.ExpVariables  = p;
%     data.trialType     = p.trialType(ii);
%     data.span          = p.spans(ii,:);
%     data.masking       = p.masking;
%     data.DAQin         = data_DAQ;
%     save([p.savedFiledir 'Trial' num2str(ii)], 'data')
%     audiowrite([p.savedWavdir 'Trial' num2str(ii) '_headOut.wav'], data.signalOut, p.postProcSRate)
%     audiowrite([p.savedWavdir 'Trial' num2str(ii) '_micIn.wav'], data.signalIn, p.postProcSRate)
%     
    color = chkRMS(data); %How loud were they?    
    set(rec, 'Color', color); set(rec, 'FaceColor', color);
    set(rec, 'Visible','on');  
    
    pause(2.0)
    set(rec, 'Visible','off');
end
close all

if bVis == 1
    OST_MULT = 500; %Scale factor for OST
    visSignals(data, 16000, OST_MULT, savedWavdir)
end
end

function s = initPerturb

s = daq.createSession('ni');
addAnalogOutputChannel(s,'Dev3',0,'Voltage'); %Output to Perturbatron
addAnalogInputChannel(s,'Dev3',0,'Voltage');  %Input from Force Sensor
addAnalogInputChannel(s,'Dev3',1,'Voltage');

s.Rate = 8000;
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

function trialType = orderTrialsLaced(numTrial, per)
%This function organizes trials in sets based on variable 'per'. 
%For example if per = 0.25, then a set of four trials will have 
%three control trials, and one catch trial. The catch trial will always be 
%the second or third trial in the set. This will generate enough sets for
%all trials in 'numTrial'.

trialType = [];
for ii = 1:(numTrial*per)
    place = round(rand) + 2; %This will break if per!= 0.25 **Fix this**
    set = zeros(1,(1/per));
    set(place) = 1;
    
    trialType = cat(2, trialType, set);
end
end

function [sigs, spans] = createPerturbs(s, numTrial, trialLen, trialType)
%This function creates the digital signal for the NIDAQ needed for each
%trial. It also keeps track of the time points when the perturbation should
%activate and deactivate. This function calculates these points for the 
%sampling rate of the NIDAQ (s.Rate). 'spans' is then converted to the 
%Audapter sampling rate (sRate) in the main function.

%The perturbation starts at a semi-random time between 1.7s and 2.1s after 
%phonation. The pertrubation stops at a semi-random time between 
%0.7 and 1.1s after the start of pertrubation.

%s:         NIDAQ object handle for keeping track of variables and I/O
%numTrial:  The number of trials
%trialLen:  The length of each trial in points
%trialType: Vector (length = numTrial) of order of trials (control/catch)

%sigs:  Per-trial digital signal to be outputted to NIDAQ
%spans: Per-trial pertrubation start and stop points to be aligned with mic data

sigs  = zeros(trialLen, numTrial);
spans = zeros(numTrial,2);
for i = 1:numTrial
    St   = round(s.Rate*(1.7 + (2.1-1.7)*rand)); %Hardset (1.7-2.1 seconds)
    pLen = round(s.Rate*(0.7 + (1.1-0.7)*rand)); %Hardset (0.7-1.1 seconds)   
    Sp   = St+pLen;    
    span = St:Sp; 

    sig  = zeros(trialLen,1);
    if trialType(i) == 1
        sig(span) = 3;
    end
    
    sigs(:,i)  = sig;
    spans(i,:) = [St Sp];
end
end

function plotPerturb(s,lenT,sig)

t = (0:1:lenT-1)/s.Rate;
plot(t, sig);

xlabel('Time'); 
ylabel('Voltage'); 
legend('Analog Output 0');
end

function [H1, H2, rec] = createVisualFB()
%Overlays for the experiment. 
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
%Feedback for the participant on how loud they are. 
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