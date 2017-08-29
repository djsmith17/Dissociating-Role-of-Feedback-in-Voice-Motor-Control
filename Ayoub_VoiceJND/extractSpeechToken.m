function baselineTokenW = extractSpeechToken(dirs)

%%making voice file for JND from prerecorded stimuli written by EHM
%%07/01/2017
%Edited 08/23/2017: Dante Smith
dirs.SavFileDir = fullfile(dirs.RecData, 'Pilot0', 'Run3', ['Pilot0' 'Run3DRF.mat']);

load(dirs.SavFileDir);
thisData  = DRF.rawData(9);      % Take the 9th trial. It will be a control trial
fs        = DRF.expParam.sRateAnal;
sample  = thisData.signalIn;   % Grab the microphone channel.
  
auto        = 1;
tokenL      = .5; %stimulus duration to be played in the JND
tokenLP     = tokenL*fs;
riseTime    = .05;
fallTime    = .05;
riseTimeP   = riseTime*fs;
fallTimeP   = fallTime*fs;
riseQperiod = (4*riseTime)^-1;
fallQperiod = (4*fallTime)^-1;
sampLen     = length(sample);
recDuration = sampLen/fs;
time        = linspace(0, recDuration, sampLen);

Words       = {'eee'};
keepaudio   = [];
Happy = 0;

if auto == 1
    stT = 2.0;
    ix1 = fs*stT;
    ix2 = ix1 + tokenLP;
else
    figure
    plot(time, sample, 'b'); ylim([-1 1])
    
    [x, y] = ginput(1);
    ix1 = round(x(1)*fs); %Choose a single point on the line with roughly .5s following it
    ix2 = ix1 + tokenLP;
end

baselineToken = sample(ix1:ix2);
window = ones(1, tokenLP);
window(1:riseTimeP) = sin(2*pi*riseQperiod*linspace(0, riseTime, riseTimeP)).^2;
window(tokenLP-fallTimeP + 1:tokenLP) = cos(2*pi*fallQperiod*linspace(0, fallTime, fallTimeP)).^2;

Window = [sin(2*pi*linspace(0,riseTime,floor(riseTime*fs))*((4*riseTime)^-1)).^2  ones(1,floor((tokenL - riseTime - fallTime) * fs)) cos(2*pi*linspace(0,fallTime,floor(fallTime*fs))*((4*fallTime)^-1)).^2];
baselineTokenW = baselineToken.*Window;                    
    
   
while (Happy == 0)
    plotH = figure;

    plot(time, sample, 'b');
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
            
            if x(2)-x(1) > tokenL  %if x(4)-x(2) > tokenL
                hold on;
                plot(linspace(t1,t2,length(sample(ix1:ix2))),sample(ix1:ix2),'r');
                Good = input('Is the selected signal good enough? (0,1):');
                if Good  == 1
                    audioSignal = sample(ix1:ix2-1+tokenL*fs); %can throw error if signal selected is too long
                    % audioSignal = sample(ix1:ix2-1+tokenL);
                    
                    %create a cosine ramp at beginning and end so it doesn't
                    %click
                    Window = [sin(2*pi*linspace(0,riseTime,floor(riseTime*fs))*((4*riseTime)^-1)).^2  ones(1,floor((tokenL - riseTime - fallTime) * fs)) cos(2*pi*linspace(0,fallTime,floor(fallTime*fs))*((4*fallTime)^-1)).^2];
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
                    plot(time,sample,'b');
                    ylim([-1 1])
                    title(sprintf('total Duration is %f',recDuration),'color',[1 1 1]);
                end
            else
                plot(time,sample,'b');
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