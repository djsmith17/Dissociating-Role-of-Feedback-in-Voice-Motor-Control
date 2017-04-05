function dfRunf0Acuity(varargin)
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

dirs = dfDirs(acuVar.project);
dirs.RecFileDir = fullfile(dirs.RecData, participant);

if exist(dirs.RecFileDir, 'dir') == 0
    mkdir(dirs.RecFileDir);
end
dirs.RecFileDir = fullfile(dirs.RecFileDir, [participant '_f0Acuity.mat']);

% Standardized values
fs           = 44100;
time         = 0:1/fs:2;
timeLen      = length(time);
freqDef      = 440;                   %Default tone frequency 69 MIDI = 440 Hz
soundwav_Def = cos(2*pi*freqDef*time);%Default tone signal 
taper        = tukeywin(timeLen, 0.05);

dist     = .4; %sets initial distortion value. This means that on trial 1, the non-reference phase will be 50%  different than the refernce phase
upstep   = .2; %initial setting that dist will go up by if user answers 'same'. These values change later based on trial number
udRatio  = 2.4483; %MacMillan 2004
downstep = upstep/udRatio; %initial setting that dist will go down by if user answers 'different'. These values change later based on trial number

MaxReversals = 18;% number of answer switches the user must make before the game ends
reversals    = 0; % counter for number of changes (reversals) over time
trial        = 0; % counter for each trial 
randdiff     = 8; % 1 in 8 trials are a catch trial

trialTypes = []; % Saves whether it is order 1 or order 2, which indicates whether the reference frequency is first or second
distVals   = []; % Saves the dist values
responses  = []; % Saves the responses
matches    = []; % Saves the user result
revValues  = []; % Saves the reversal values, will quit at 12

changeDirection = 1; %Placeholder to remember "direction" of delta changes to count reversals
correctInARow   = 0; %Tracks # of consecutively correct trials (regardless of catch or different)
                     %Resets after 2 correct in a row

s = RandStream.create('mt19937ar', 'seed', sum(100*clock)); %make it so randperm doesn't call the same thing everytime when matlab is opened
RandStream.setGlobalStream(s);

%Start her up
disp('Hit return when ready to start')
pause

% Loop until we get the necessary # of reversals
while reversals < MaxReversals 
    trial = trial + 1; % advance the trial every time we go through this loop
    
    RandomCatchTrial = randperm(randdiff, 1); % 1/8 if the time it is the same and 50% it is different
    
    % Create the comparison stimulus
    if RandomCatchTrial < 8
        diffTokens = 1; 
        trialType  = randperm(2,1); %different trials, order 1 reference is first, order 2 reference is second 
        distVal    = dist;
    else
        diffTokens = 0;
        trialType  = 3; 
        distVal    = 0;
    end
    trialTypes = cat(1, trialTypes, trialType);
    distVals   = cat(1, distVals, distVal);
    
    newf = freqDef*2^(dist/12);
    if trialType == 1
        token1 = soundwav_Def;
        token2 = cos(time*2*pi*newf);
    elseif trialType == 2
        token1 = cos(time*2*pi*newf);
        token2 = soundwav_Def;
    elseif trialType == 3
        token1 = soundwav_Def;
        token2 = soundwav_Def;
    end
    
    token1_proc = preProcessToken(token1, taper);
    token2_proc = preProcessToken(token2, taper);
    
    %The actual trial. Listen to two tokens. Are they the same or different
    disp('Hit return to begin next trial')% what the users sees on the screen
    pause
    sound(token1_proc, fs) %Token 1
    pause(2)
    sound(token2_proc, fs) %Token 2
    fprintf('Same (0) or Different (1)?');

    YourAnswer = getkey; %getkey is a function that finds the value of the next keypress, getkey MUST be in the same folder as this code
    fprintf(' %c\n', YourAnswer) %displays the key that was entered

    if YourAnswer == 49 % key press = 1
        response = 1; %the response was 1 (different)
    elseif YourAnswer == 48 % key press = 0
        response = 0; %the response was 0 (same)
    end
    responses = cat(1, responses, response);
    
    match = isequal(diffTokens, response); %Was the subject correct?
    matches = cat(1, matches, match);
    
    for Phase = 1:2 %phase 1 plays first then phase 2
        if (trialType == 3) || (trialType == 1 && Phase == 1) || (trialType == 2 && Phase == 2)
           %if catchtrial, both phases play reference
           %OR if different trial order 1 AND first phase play reference
           %OR if different trial order 2 AND second phase play reference
            
           if trialType == 1 
               distVals = [distVals dist]; %saves all the dist values BEFORE the person hears them/makes a decisions
           end

           soundwav = soundwav_Def; %plays the baseline - doesn't change
        else % makes the comparision stimuli
            if trialType == 2
                distVals = [distVals dist];
            end   %saves all the dist values BEFORE the person hears them/makes a decisions

            newf = freqDef*2^(dist/12);
            soundwav = cos(time*2*pi*newf);
        end

       % Window the beginning and end of the stimuli slightly to avoid clicks
        soundwav = soundwav .* taper';

        if Phase == 1
            disp(['Hit return to begin next trial'])% what the users sees on the screen
            pause
        end

        sound(.05*soundwav*10^(-6/20), fs); %65 dB on MOTU equipment

        if Phase == 1
            pause(2) % Play the two sounds in succession (with a short pause between)
        end

        if Phase == 2
            fprintf('Same (0) or Different (1)?'); % ask them if they were the same or different
            YourAnswer = getkey; %getkey is a function that finds the value of the next keypress, getkey MUST be in the same folder as this code
            fprintf(' %c\n', YourAnswer) %displays the key that was entered
            
            if YourAnswer == 49 % key press = 1
                response(trial) = 1; %the response was 1 (different)
            elseif YourAnswer == 48 % key press = 0
                response(trial) = 0; %the response was 0 (same)
            end

            % Determine if the user was correct or incorrect

            % If it is a catch trial, then isDiff will be 0, and 
            % a response of "Same" (0) will be a match. 

            % If it is not a catch trial, isDiff will be 1, and
            % a response of "Different" (1) will be a match.

            % In all other cases the trial will not be a match
            match(trial) = isequal(diffTokens,(response(trial))); 
            
            if dist < 0
              dist     = 0.02;
              upstep   = 0.01;
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
            if (trial > 1) %nothing happens on trial 1, it does not get included into 2 up and 1 down calculations
                if (match(trial)) %this is if match = 1, so they were correct this trial
                    if (correctInARow == 1) %this means they were correct last trial as well
                        dist = dist - downstep; %they now have got 2 correct in a row, so move the dist down to make it closer to the reference, aka harder

                        if (changeDirection == 1)% if it did change direction count a reversal, record the dist value and then note that it already changed directions
                            reversals = reversals + 1; % counting the reversals
                            revValues = [revValues dist]; %recording the revValues of the dist that JUST happened 
                            changeDirection = -1; %nothing that it just changed directions
                        end
                        correctInARow = 0; % reset the # in a row counter after you have changed the dist
                    else
                        dist = dist; % no change in dist
                        % This is the first correct in a row.  No change in dist
                        correctInARow = 1; % mark our # in a row counter
                    end
                else
                    % The subject got the last trial wrong -> increase delta
                    dist = dist + upstep;
                    if (changeDirection == -1)
                        % This was a reversal in direction of changes
                        reversals = reversals + 1;
                        revValues = [revValues dist];
                        changeDirection = 1;
                    end
                    correctInARow = 0;
                end
            else
                % This was the first trial of the experiments - do nothing
                dist(trial) = dist(trial);
            end

            if dist < 0
                dist = .02; %disp('it was negative!')
            end
            if trial > 84
                reversals = MaxReversals + 1;
            end
        
            ResultMatrix = [distVals', matches', responses', trialTypes'];
            save(dirs.SavFileDir, 'ResultMatrix', 'revValues');
        end
    end
end

Last5Mean = mean(revValues(end-5:end));
end

function token_proc = preProcessToken(token, taper)
A  = .05;
xp = -6/20; 

token_tape = token .* taper';
token_proc = A*token_tape*10^(xp);
end

% function drawJNDResults()
% figure;
% plot(delta,'-o');
% xlabel('Trial number');
% ylabel('\Delta f (Hz)','Interpreter','Tex');
% 
% Calculate the average delta over the last nJNDTrials trials
% JND = mean(delta(end-nJNDTrials:end));
% msgbox(sprintf('JND estimate is %2.2f Hz',JND));
% end