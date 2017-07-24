function meanJNDVal = dfRunf0Acuity(varargin)
%Edit 08-15-2013: Cara Stepp, amplitude changes
%Edit 06-13-2014: Liz Heller Murray
%Edit 04-07-2015: Defne Abur, normalized sound, 65dB on MOTU
%Edit 04-09-2015: Defne Abut, added 'order', 'ordersave' variables
%Edit 05-21-2015: Liz Heller Murray, only change when two correct in a row
%Edit 02-26-2017: Dante Smith, just some general organization

%% GAME INITIALIZATION
if isempty(varargin)
    participant = 'null'; 
else
    participant = varargin{1};
end

acuVar.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
acuVar.participant  = participant;
acuVar.run          = 'Run1';

%dfDirs is a separate function that regulates my dirs and paths. 
%dirs is a structure with my relevant paths for data files. 
dirs = dfDirs(acuVar.project);
dirs.RecFileDir = fullfile(dirs.RecData, participant);
dirs.SavFileDir = fullfile(dirs.RecData, 'Pilot0', 'Run3', ['Pilot0' 'Run3DRF.mat']);

if exist(dirs.RecFileDir, 'dir') == 0
    mkdir(dirs.RecFileDir);
end
dirs.RecFileDir = fullfile(dirs.RecFileDir, [participant '_f0Acuity' acuVar.run '.mat']);

%This should return DRF
load(dirs.SavFileDir);
thisData  = DRF.rawData(8);      % Take the 8th trial. It will be a control trial
speech    = thisData.signalIn;   % Grab the microphone channel.
baseToken = speech(16000:47999); % 2s long sample


levels = (1:10)*0.02;
modTokens = generateSpeechTokens(dirs, baseToken, levels);

% Standardized values
fs           = 44100;
time         = 0:1/fs:2;
timeLen      = length(time);
freqDef      = 440;                     %Default tone frequency 69 MIDI = 440 Hz
soundwav_Def = cos(2*pi*freqDef*time);  %Default tone signal 
taper        = tukeywin(timeLen, 0.05); %Windowing function

dist     = .4; %sets initial distortion value. This means that on trial 1, the non-reference phase will be 50%  different than the refernce phase
upstep   = .2; %initial setting that dist will go up by if user answers 'same'. These values change later based on trial number
udRatio  = 2.4483; %MacMillan 2004
downstep = upstep/udRatio; %initial setting that dist will go down by if user answers 'different'. These values change later based on trial number

MaxReversals = 5;% number of answer switches the user must make before the game ends
reversals    = 0; % counter for number of changes (reversals) over time
trial        = 0; % counter for each trial 
randdiff     = 8; % 1 in 8 trials are a catch trial

trialTypes = []; % Saves whether it is order 1 or order 2, which indicates whether the reference frequency is first or second
distVals   = []; % Saves the dist values
responses  = []; % Saves the responses
matches    = []; % Saves the user answer
revValues  = []; % Saves the reversal values, will quit at 12

changeDirection = 1; %Placeholder to remember "direction" of delta changes to count reversals
correctInARow   = 0; %Tracks # of consecutively correct trials (regardless of catch or different)
                     %Resets after 2 correct in a row

%Make it so randperm doesn't call the same thing everytime when matlab is opened
s = RandStream.create('mt19937ar', 'seed', sum(100*clock)); 
RandStream.setGlobalStream(s);

%Start her up
disp('Hit return when ready to start')
pause

% Loop until we get the necessary # of reversals
while reversals < MaxReversals 
    trial = trial + 1; % advance the trial every time we go through this loop 
    RandomCatchTrial = randperm(randdiff, 1); % 1/8 if the time it is the same and 50% it is different
    
    %Determine trial type and record dist to be used
    if RandomCatchTrial < randdiff
        diffTokens = 1; 
        trialType  = randperm(2,1); %different trials, order 1 reference is first, order 2 reference is second 
        distVal    = dist;
    else
        diffTokens = 0;
        trialType  = 3; 
        distVal    = 0;
    end
    trialTypes = cat(1, trialTypes, trialType); %Save trialTypes
    distVals   = cat(1, distVals, distVal);     %Save the distVal
    
    %Set up the audio tokens to be played
    freqNew = freqDef*2^(dist/12);
    if trialType == 1
        token1 = soundwav_Def;
        token2 = cos(2*pi*freqNew*time);
    elseif trialType == 2
        token1 = cos(2*pi*freqNew*time);
        token2 = soundwav_Def;
    elseif trialType == 3
        token1 = soundwav_Def;
        token2 = soundwav_Def;
    end
    
    %Process the audio tokens for easy perception. Func is below
    token1_proc = preProcessToken(token1, taper); 
    token2_proc = preProcessToken(token2, taper);
    
    %The actual trial. Listen to two tokens. Are they the same or different
    fprintf('Hit return to begin next trial\n')% what the users sees on the screen
    pause
    sound(token1_proc, fs) %Token 1
    pause(2)
    sound(token2_proc, fs) %Token 2
    
    %Ask the subject if the tones were the same or different
    fprintf('Same (0) or Different (1)?');    
    YourAnswer = getkey; %getkey is a function that finds the value of the next keypress, getkey MUST be in the same folder as this code
    fprintf(' %c\n', YourAnswer) %displays the key that was entered

    if YourAnswer == 49 % key press = 1
        response = 1; %the response was 1 (different)
    elseif YourAnswer == 48 % key press = 0
        response = 0; %the response was 0 (same)
    end
    responses = cat(1, responses, response); %Save response results
    
    %Determine if the subject was correct
    match = isequal(diffTokens, response); %Was the subject correct?
    matches = cat(1, matches, match);      %Save match results

    if dist < 0
        dist   = 0.02;
        upstep = 0.01;
    else
        if trial > 10 %reduce step size
            upstep = 0.1;
        elseif trial > 20 %reduce step size again
            upstep = 0.05;
        elseif trial > 30 %reduce step size again (last time)
            upstep = 0.02;
        end
    end
            
    downstep = upstep/udRatio;

    % Adaptively adjust dist as required by 2-down, 1-up
    if trial > 1 %Trial 1 is not included in the 2-down and 1-up calculations
        if match == 1 %The subject was correct -> maybe decrease dist
            if correctInARow == 1 %Subject was correct last trial as well, time to move dist
                correctInARow = 0;      %Reset correctness counter
                dist = dist - downstep; %Decrease dist (Task is harder)

                if changeDirection == 1 %Reversal flipping to decreasing dist
                    reversals = reversals + 1;           %Record # of reversals
                    revValues = cat(1, revValues, dist); %Record dist at the reversal
                    changeDirection = -1;                %Use decreasing values of dist until subject answers incorrectly
                end               
            else
                correctInARow = 1; %First correct in a row.  No change in dist                
            end
        elseif match == 0 %The subject was wrong -> absolutely increase dist
            correctInARow = 0;    %Reset correctness counter
            dist = dist + upstep; %Increase dist (Task is easier)
            
            if changeDirection == -1 %Reversal flipping to increasing dist
                reversals = reversals + 1;           %Record # of reversals
                revValues = cat(1, revValues, dist); %Record dist at the reversal
                changeDirection = 1;                 %Use increasing values of dist until subject answers correctly
            end         
        end
    end

    if dist < 0
        dist = .02; disp('it was negative!')
    end
    if trial > 84
        reversals = MaxReversals + 1;
    end
end

f0AcuityResults = [trialTypes, distVals, responses, matches];
save(dirs.RecFileDir, 'f0AcuityResults', 'revValues');

meanJNDVal = mean(revValues);

drawJNDResults(dirs, acuVar, f0AcuityResults, revValues, meanJNDVal)
end

function token_proc = preProcessToken(token, taper)
A  = .05;
xp = -6/20; 

token_tape = token .* taper';
token_proc = A*token_tape*10^(xp);
end

function modTokens = generateSpeechTokens(dirs, baseToken, levels)
modTokens = [];

for n = 1:length(levels)
    modifiedToken = dfGenerateOfflineTokens(baseToken, dirs, 'male', levels(n));
    modTokens = cat(1, modTokens, modifiedToken');
end
end

function drawJNDResults(dirs, acuVar, f0AcuityResults, revValues, meanJNDVal)

ind = find(f0AcuityResults(:,2) ~= 0);
butts = f0AcuityResults(ind,:);
numTrial = length(butts);
trialVec = 1:numTrial;
distVec  = butts(:,2);

f0AcuityFig = figure('Color', [1 1 1]);
plot(trialVec, distVec)
xlabel('Trial #', 'FontSize', 18, 'FontWeight', 'bold')
ylabel('Distortion Val', 'FontSize', 18, 'FontWeight', 'bold')
title(['f0 Acuity Results for ' acuVar.participant ' ' acuVar.run], 'FontSize', 18, 'FontWeight', 'bold')
box off

end