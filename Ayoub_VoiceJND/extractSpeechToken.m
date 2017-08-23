function baselineToken = extractSpeechToken(dirs)

%%making voice file for JND from prerecorded stimuli written by EHM
%%07/01/2017
%Edited 08/23/2017: Dante Smith
clear all
close all
project = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
dirs = dfDirs(project);
dirs.SavFileDir = fullfile(dirs.RecData, 'Pilot0', 'Run3', ['Pilot0' 'Run3DRF.mat']);

load(dirs.SavFileDir);
thisData  = DRF.rawData(9);      % Take the 9th trial. It will be a control trial
fs        = DRF.expParam.sRateAnal;
Y_sample  = thisData.signalIn;   % Grab the microphone channel.

% [Y_sample, fs] = audioread('TestVoice.wav');
  
StimulusDur = .5; %this is the duration of the stimulus that will be played in the JND
riseTime    = .05;
fallTime    = .05;
recDuration = length(Y_sample)/fs;
Words       = {'aaa'};
keepaudio   = [];

figure1 = figure('Color',[1 1 1],'Menubar','none');

drawnow;
Happy = 0;
   
while (Happy == 0)
    plotH = figure;
    times = linspace(0, recDuration, length(Y_sample));
    times = times';
    plot(times, Y_sample, 'b');
    ylim([-1 1])
    % title(sprintf('total Duration is %f',Record_Time));
    Happy = input('Are you happy with the signal (0,1):');
    if Happy == 1
        
        Select = 0;
        while Select == 0
            [x, y] = ginput(2);
            ix1 = round(x(1)*fs);   %ix1 = round(x(2)*44100);
            ix2 = round(x(2)*fs);   %ix2 = round(x(4)*44100);
            t1 = ix1/fs;
            t2 = ix2/fs;
            
            if x(2)-x(1) > StimulusDur  %if x(4)-x(2) > StimulusDur
                hold on;
                plot(linspace(t1,t2,length(Y_sample(ix1:ix2))),Y_sample(ix1:ix2),'r');
                Good = input('Is the selected signal good enough? (0,1):');
                if Good  ==1
                    audioSignal = Y_sample(ix1:ix2-1+StimulusDur*fs); %can throw error if signal selected is too long
                    % audioSignal = Y_sample(ix1:ix2-1+StimulusDur);
                    
                    %create a cosine ramp at beginning and end so it doesn't
                    %click
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
                    %     audioSignal = .5 * (audioSignal ./ (rms(audioSignal))); %% was used to scale
                    %     amplitude up. Currently results in clippling so
                    %     need to test if want to add back in
                    audiowrite('voice_use.wav',audioSignal,44100);
                    
                    Select = 1;
                    clf(plotH);
                    close (plotH);                   
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
        hold off       
    end 
end