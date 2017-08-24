function out = GAP_PitchJND_Passive_500()

close all;
ET = tic;
rng('shuffle');
%% define experimental parameters
%edited 10/20/2016 DA

calib_check = input('Has this microphone position been calibrated? (0 = no,1 = yes):');

if calib_check ==0
    errordlg('Run the IntensityCalibration script');
    return
end

cue_check = input('Did you check the CueMix Fx configuration is correct? (0 = no,1 = yes):');
if cue_check == 0
    errordlg('Check CueMix setting!');
    return
end

prompt={'CueMix Mic Trim: '}; %ask for cuemix trim settings. Put Trim number seen in Mic inputs. 
name='CueMix';
numlines=1;
defaultanswer={'00'}; %input box for cuemix trim setting
answer=inputdlg(prompt,name,numlines,defaultanswer);

cuemixMic = str2num(answer{1}); %cue mix trim trim setting

prompt={'CueMix MainOut (100 for max, otherwise clock orientation): '}; %This line asks for cuemix trim settings for MainOut. MainOut doesn't have values, so max settings is "100", otherwise, specify clock orientation (for example, for the setting level with the "12" on the bar to the left, it would be 9). 
name='MainOut';
numlines=1;
defaultanswer={'00'}; %input box for cuemix trim setting
answer=inputdlg(prompt,name,numlines,defaultanswer);
cuemixMainout = str2num(answer{1});%main out setting 



prompt={'Subject ID:',...
    'Session ID:',...
    'Group:',...
    'Gender ("male" or "female")'};
name='Subject Information';
numlines=1;
defaultanswer={'CTGAP','Pitch_JND_Passive_100ms','Control','female'};
answer=inputdlg(prompt,name,numlines,defaultanswer);

if isempty( answer )
    return
end

num_trials = questdlg('Practice or Full?','Length','Practice','Full','Full') ;
switch num_trials
    case 'Practice'
        num_trials = 'Practice';
        UD.totalTrials = 9;
    case 'Full'
        num_trials = 'Full';
        UD.totalTrials = 60; %max number of trials if max trials/reversals not reached
end

subjectID = answer{1};
session = answer{2};
group =  answer{3};
Gender =  answer{4};
% create the folder path to save the data files
baseFilename = ['data\' group, '\', subjectID, '\', session, '\',num_trials,'\'];

% check if the foler exists (to avoid overwrite)
if ~exist( baseFilename , 'dir' )
    mkdir( baseFilename ) ;
    while exist(baseFilename,'dir') ~= 7
            mkdir(baseFilename);
    end
else
    overwrite = inputdlg({'File already exists! Do you want to overwrite?'},'Overwrite',1,{'no'});
    if ~strcmp(overwrite,'yes') & ~strcmp(overwrite,'YES') & ~strcmp(overwrite,'Yes')
        return;
    end
end

%% recording audio samples
Words = {'aaa'};
recDuration = 4; %this is the duration of the audio recording of /a/ production
StimulusDur = .5; %this is the duration of the stimulus that will be played in the JND
riseTime = .05; %cosine square (ramp) to avoid click
fallTime = .05; %cosine square (ramp) to avoid click
initializeAudapter (Gender); %this does not change anything for pitch but it does for formant

Audapter('reset');
Audapter('start');
pause(1);
Audapter('stop');

figure1 = figure('Color',[0 0 0],'Menubar','none');
textBox1 = annotation(figure1,'textbox',...
    [0.25 0.65 0.5 0.25],...
    'Color',[1 1 1],...
    'String','READY',...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'FontSize',130,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[0 0 0],...
    'Visible','on');
drawnow;
Happy =0;
pause(5); %to give time to move window to other screen
% Present stimuli and collect the results
while (Happy ==0)
    clf(figure1);
    textBox1 = annotation(figure1,'textbox',...
        [0.25 0.65 0.5 0.25],...
        'Color',[1 1 1],...
        'String','READY',...
        'LineStyle','none',...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle',...
        'FontSize',130,...
        'FontName','Arial',...
        'FitBoxToText','off',...
        'EdgeColor','none',...
        'BackgroundColor',[0 0 0],...
        'Visible','on');
    set(textBox1,'String',Words{1})
    Audapter('reset');
    Audapter('start');
    pause(recDuration);
    Audapter('stop');
    clf(figure1);
    % save the data
    data1 = AudapterIO('getData');
    plotH = figure;
    Y_sample = data1.signalIn;
    Fs = 16000;
    times = linspace(0,recDuration,length(data1.signalIn));
    plot(times,Y_sample,'b');
    ylim([-1 1])
    title(sprintf('total Duration is %f',recDuration));
    Happy = input('Are you happy with the signal (0,1):');
    if Happy == 1
        
        Select = 0;
        while Select ==0
            [x,y] = ginput(4);
            ix1 = round(x(2)*16000);
            ix2 = round(x(4)*16000);
            t1 = ix1/16000;
            t2 = ix2/16000;
            
            if x(4)-x(2) > StimulusDur  %* Fs
                hold on;
                plot(linspace(t1,t2,length(Y_sample(ix1:ix2))),Y_sample(ix1:ix2),'r');
                Good = input('Is the selected signal good enough? (0,1):');
                if Good  ==1
                    audioSignal = Y_sample(ix1:ix2-1+StimulusDur*Fs);
                    %                     Window = [sin(2*pi*linspace(0,riseTime,round(riseTime*Fs))*((4*riseTime)^-1)).^2  ones(1,round((StimulusDur - riseTime - fallTime) * Fs)) cos(2*pi*linspace(0,fallTime,round(fallTime*Fs))*((4*fallTime)^-1)).^2];
                    audioSignal = audioSignal(:);%.*Window(:);
                    Select = 1;
                    close (plotH);
                    %                     plot(audioSignal)
                else
                    Good = 0;
                    Select = 0;
                    plot(times,Y_sample,'b');
                    ylim([-1 1])
                    title(sprintf('total Duration is %f',recDuration),'color',[1 1 1]);
                end
            else
                plot(times,Y_sample,'b');
                ylim([-1 1])
                title(sprintf('total Duration is %f',recDuration),'color',[1 1 1]);
                hhError = errordlg('the selected signal is too short!');
                %                 pause(2)
                %                 close(hhError)
            end
        end
    else
        close (plotH)
    end
end
close all;
% Normfact = sqrt(mean(audioSignal.^2));%ADDED BY CARA TO TRY TO NORMALIZE SOUND AMPLITUDES!
% audioSignal = .5* (audioSignal./Normfact);

% audioSignal = .5 * (audioSignal ./ (max(abs(audioSignal))));
audioSignal = .5 * (audioSignal ./ (rms(audioSignal)));


%% Setting up the up-down paradigm (modified based on Palam)

UD.up = 1;    % Number of consecutive responses before an increase
UD.down = 2;  % Number of consecutive responses before a decrease
stepSize = 4; %This is something to tune; in cents
UD.stepSizeUp = stepSize / 1; %Levitt (1971) 2/1 rule for 71% in MacMillian Chapter 11 with step per Garcia-Perez (1998); Was: Size of step up ; stepSize/ .5488 ensures 80.35 % correct; see Garcia-Perez 1998
UD.stepSizeDown = stepSize; % Size of step down
UD.stopCriterion = 'reversals'; % stop the procedure based on number of 'trials' | 'reversals'
UD.stopRule = 10;  %stop procedure after this number of trials/reversals
UD.startValue = 50; % initial difference in cents between speaker's fo and fo of stimulus in headphones
UD.xMax = 200; %max difference between speaker's fo and fo of stimulus in headphones
UD.xMin = 0; %min difference between speaker's fo and fo of stimulus in headphones
UD.truncate = 'yes';
UD.response = [];
UD.stop = 0;
UD.u = 0;
UD.d = 0;
UD.direction = [];
UD.reversal = 0;
UD.xCurrent = UD.startValue;
UD.x = UD.startValue;
UD.xStaircase = [];
% UD.totalTrials = 60; %max number of trials if max trials/reversals not reached
waitForKeyPress = 3 ; % in seconds
baseTime = .5; %this is the interstimulus interval within each trial (between stimulus 1 and stimulus 2 in pair) in seconds
%% visual presentation
monitorSize = get(0,'Monitor');
if size(monitorSize,1) == 1
    figPosition = [1 200 monitorSize(3) monitorSize(4)-200];
elseif size(monitorSize,1) == 2
    figPosition = [monitorSize(2,1) monitorSize(2,2)+20 monitorSize(2,3) monitorSize(2,4)];
end

figure1 = figure('Color',[0 0 0],'Menubar','none','Position',[-1280 1050 1281 1026]);

h2 = annotation(figure1,'textbox',...
    [0.38 0.46 0.2 0.2],...
    'Color',[1 1 1],...
    'String','READY',...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'FontSize',130,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[0 0 0],...
    'Visible','on');

h3 = annotation(figure1,'textbox',...
    [0.025 0.15 0.45 0.3],...
    'String',{'< DIFFERENT'},... %was 'YES'
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'FontSize',60,...
    'FontName','Arial',...
    'LineStyle','none',...
    'BackgroundColor',[1 1 1],...
    'Color',[0 0 0],...
    'Visible','off');

h4 = annotation(figure1,'textbox',...
    [0.52 0.15 0.45 0.3],...
    'String',{'SAME >'},... %was 'NO'
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'FontSize',60,...
    'FontName','Arial',...
    'LineStyle','none',...
    'BackgroundColor',[1 1 1],...
    'Color',[0 0 0],...
    'Visible','off');

drawnow;
pause(5);

tr = 0;
while (UD.stop == 0) & tr < UD.totalTrials
    tr = tr +1;
    %Present the word
    set(h2,'String','+')
    drawnow;
    
    tempVar = randperm(5);
    Pert = UD.xCurrent/100;
    if tempVar(1) == 1 || tempVar(1) == 3   % scenario I (first one is Pert) : % 40% of trials
        offline_pitch_JND_passive (audioSignal,Pert,Gender,StimulusDur,riseTime,fallTime,baseTime);
        %         pause(baseTime);
        offline_pitch_JND_passive (audioSignal,0.01,Gender,StimulusDur,riseTime,fallTime,0); %.01 perturbation has been applied to reference stimulus to eliminate differences in signal quality related to Audapter processing
        conVar = 1;
    elseif tempVar(1) == 2 || tempVar(1) == 4 % scenario II (first one is no Pert) : % 40% of trials
        offline_pitch_JND_passive (audioSignal,0.01,Gender,StimulusDur,riseTime,fallTime,baseTime); %.01 perturbation has been applied to reference stimulus to eliminate differences in signal quality related to Audapter processing
        %         pause(baseTime);
        offline_pitch_JND_passive (audioSignal,Pert,Gender,StimulusDur,riseTime,fallTime,0);
        conVar = 1;
    else % Catch trials : 20% of trials
        offline_pitch_JND_passive (audioSignal,0.01,Gender,StimulusDur,riseTime,fallTime,baseTime); %.01 perturbation has been applied to reference stimulus to eliminate differences in signal quality related to Audapter processing
        %         pause(baseTime);
        offline_pitch_JND_passive (audioSignal,0.01,Gender,StimulusDur,riseTime,fallTime,0); %.01 perturbation has been applied to reference stimulus to eliminate differences in signal quality related to Audapter processing
        conVar = 0; %catch trials will be randomly presented but will not be included in the adaptive procedure
        
    end
    
    %% Present the YES/NO question
    %     set(h1,'Visible','off');
    set(h2,'String','PITCH','FontSize',80)%was 'WERE THEY DIFFERENT'
    set(h3,    'Visible','on');
    set(h4,    'Visible','on');
    drawnow
    keyCorrect = 1;
    while keyCorrect
        
        keyCorrect = 0;
        [bb, ReactionTime(tr)] = GetKey_Ayoub(1,waitForKeyPress);
        
        % wait until a correct key is pressed
        if (bb ~= 28) & (bb ~= 29)
            keyCorrect = 1;
        end
        if isempty(bb)
            keyCorrect = 1;
        end
        
        if bb == 28            %28 is DIFFERENT" | Left ARROW KEY ; was YES
            response = 1;
        elseif bb == 29        %29 is "SAME"  | Right ARROW KEY ; was NO
            response = 0;
        end
        
    end
    set(h2,'String','','FontSize',120)
    set(h3,    'Visible','off');
    set(h4,    'Visible','off');
    drawnow
    if conVar == 1 % when it is a real trial update the UD structure
        UD = adaptiveUD_Update(UD, response);
        UD.catchResponse(tr,1) = NaN;
    else % when it is a catch trial do not update UD structure (i.e., do not change the up-down steps based on catch trials)
        UD.catchResponse(tr,1) = response;
    end
    pause(1) %this is between two trials
    
    
end
close all;
UD .audioSignal = audioSignal;
%CES - To disable saving of reactiontimes (we are doing this to avoid confusion since the investigator is keying in the selection at present)
UD.reactionTime = ones(size(ReactionTime))*10000;

FileName= [baseFilename 'responseResults.mat'];
dataFileName = [baseFilename 'dataFile.mat'];
save(FileName,'UD');

CueMixSave = [baseFilename 'CueMixSettings.mat']; %makes a new .mat file to save 
CueMixSettings.cuemixMic = cuemixMic; %defines saving of user input mic trim
CueMixSettings.cuemixMainout = cuemixMainout; %defines saving of user input main out value
save(CueMixSave,'CueMixSettings');%saves data

clc;

% thresholdAnalyzeUD(UD, 'reversals',4);
elapsed_time = toc(ET)/60;
disp (sprintf('Total time: %f (min)',elapsed_time));
switch num_trials
    case 'Practice'
        out = [];
    case 'Full'
        out = thresholdAnalyzeUD(UD, 'reversals',4)*.01;
        dataFile.time = elapsed_time;
        dataFile.score = out;
        dataFile.totalTrials = length(UD.catchResponse);
        dataFile.JNDTrials = length(UD.reversal);
        dataFile.catchTrials = length(UD.catchResponse) - length(UD.reversal);
        dataFile.reversals = max(UD.reversal);
        dataFile.catchCorrect = sum(UD.catchResponse == 0);
        save(dataFileName,'dataFile');
end
end



function data1 = offline_pitch_JND_passive (dataFile,Pert,genderSubject,StimulusDur,riseTime,fallTime,basetime)
%% CONFIG
addpath c:/speechres/commonmcode;
 cds('audapter_matlab');
audioInterfaceName = 'ASIO4ALL';%'MOTU MicroBook';%

sRate = 48000;  % Hardware sampling rate (before downsampling)
downFact = 3;
frameLen = 96;  % Before downsampling
noiseWavFN = 'mtbabble48k.wav';
% Audapter('deviceName', audioInterfaceName);
Audapter('setParam', 'downFact', downFact, 0);
Audapter('setParam', 'sRate', sRate / downFact, 0);
Audapter('setParam', 'frameLen', frameLen / downFact, 0);

%%
ostFN = '../example_data/offline_pitch_JND_passive.ost';
pcfFN = '../example_data/offline_pitch_JND_passive.pcf';

write2pcf(pcfFN , Pert)

check_file(ostFN);
check_file(pcfFN);
Audapter('ost', ostFN, 0);
Audapter('pcf', pcfFN, 0);

%%
params = getAudapterDefaultParams(lower(genderSubject));

params.f1Min = 0;
params.f2Max = 5000;
params.f2Min = 0;
params.f2Max = 5000;
params.pertF2 = linspace(0, 5000, 257);
params.pertAmp = 0.0 * ones(1, 257);
params.pertPhi = 0.0 * pi * ones(1, 257);
params.bTrack = 1;
params.bShift = 1;
params.bRatioShift = 1;
params.bMelShift = 0;
maxPBSize = Audapter('getMaxPBLen');
check_file(noiseWavFN);
[w, fs] = audioread(noiseWavFN);
if fs ~= params.sr * params.downFact
    w = resample(w, params.sr * params.downFact, fs);
end
if length(w) > maxPBSize
    w = w(1 : maxPBSize);
end
% Audapter('setParam', 'datapb', w, 1);
params.fb3Gain = 0.1;
params.fb = 1;
params.pertAmp = Pert * ones(1, 257);
params.pertPhi = 0.0 * pi * ones(1, 257);
params.fb = 1;
% AudapterIO('init', params);
params.bDetect = 1;
params.rmsThresh = 0.01;
params.bRatioShift = 1;

params.bBypassFmt = 1; %formant tracking bypassed, flag that not shifting formants during this script
params.bTrack =0; %not necesssary because formants are not being shifted in this script
%params.bShift = 1;
%params.bBypassFmt = 0;

params.bPitchShift = 1;


fs = 16000;

sigIn = dataFile;
sigIn = resample(sigIn, fs * downFact, fs);
sigInCell = makecell(sigIn, frameLen);

params.rmsClipThresh=0.01;
params.bRMSClip=1;
AudapterIO('init', params);
Audapter('setParam', 'rmsthr', 5e-3, 0);
Audapter('reset');
for n = 1 : length(sigInCell)
    Audapter('runFrame', sigInCell{n});
end
data1 = AudapterIO('getData');

%% normalize
audioSignal = data1.signalOut;
RMS_Source = sqrt(mean(dataFile.^2));
RMS_Target = sqrt(mean(audioSignal.^2));
audioSignal = audioSignal * (RMS_Source/RMS_Target);%this is normalizing the amplitude
%Window = [sin(2*pi*linspace(0,riseTime,riseTime*fs)*((4*riseTime)^-1)).^2  ones(1,(StimulusDur - riseTime - fallTime) * fs) cos(2*pi*linspace(0,fallTime,fallTime*fs)*((4*fallTime)^-1)).^2];
Window = [sin(2*pi*linspace(0,riseTime,floor(riseTime*fs))*((4*riseTime)^-1)).^2  ones(1,floor((StimulusDur - riseTime - fallTime) * fs)) cos(2*pi*linspace(0,fallTime,floor(fallTime*fs))*((4*fallTime)^-1)).^2];
Window = Window(:);%this is a cosine square
audioSignal = audioSignal(:);
if length(Window) > length(audioSignal)
    audioSignal = audioSignal(:) .* Window(1:length(audioSignal));
elseif  length(Window) < length(audioSignal)
    audioSignal = audioSignal(1:length(Window)) .* Window(:);
else
    audioSignal = audioSignal(:) .* Window(:);
end
audioSignal = audioSignal(:);


audioSignal = audioSignal(:).*Window(:);

audioSignal = .5*audioSignal/rms(audioSignal);
% playerObj = audioplayer(audioSignal, fs);
% playblocking(playerObj);
sound(audioSignal, fs)
pause(basetime+length(audioSignal)/fs+.01)

end

