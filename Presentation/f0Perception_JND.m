function f0Perception_JND()
%Edit 08-15-2013: Cara Stepp, amplitude changes
%Edit 06-13-2014: Liz Heller Murray
%Edit 04-07-2015: Defne Abur, normalized sound, 65dB on MOTU
%Edit 04-09-2015: Defne Abut, added 'order', 'ordersave' variables
%Edit 05-21-2015: Liz Heller Murray, only change when two correct in a row
%Edit 02-26-2017: Dante Smith, just some general organization

%% GAME INITIALIZATION
clear all; close all; % empty the MATLAB workspace

%generates the GUI that pops up at the beginning of the game
prompt        = {'Subject:'};
name          = 'Hearing JND';
numlines      = 1;
defaultanswer = {'s001'};
 
answer = inputdlg(prompt,name,numlines,defaultanswer);

trial = ['t001']; %will hold all of the summary information you need to analyize
outfile = fullfile(answer{1}, trial) ;%creates file 't001'

s = RandStream.create('mt19937ar','seed',sum(100*clock)); %make it so randperm doesn't call the same thing everytime when matlab is opened
RandStream.setGlobalStream(s);

if isempty( answer )
    return
end

if ~exist( answer{ 1 } , 'dir' )
    mkdir( answer{ 1 } ) ;
else
    overwrite = inputdlg({'File Exists: Are you sure you want to overwrite?'},'',1,{'no'}); %makes sure you don't write over a file you already have done
    if ~strcmp(overwrite,'yes') & ~strcmp(overwrite,'YES') & ~strcmp(overwrite,'Yes') 
      return;
    end
end


%% set intial variables

% Standardized values
fs       = 44100;
timeT    = [0:1/fs:2];
soundwav = cos(timeT*2*pi*440); %69 MIDI = 440 Hz

dist     = .4; %sets initial distortion value. This means that on trial 1, the non-reference phase will be 50%  different than the refernce phase
upstep   = .2; %initial setting that dist will go up by if user answers 'same'. These values change later based on trial number
downstep = upstep / 2.4483; %initial setting that dist will go down by if user answers 'different'. These values change later based on trial number


MaxReversals  = 18;  %number of answer switches the user must make before the game ends
correctInARow = 0; % keep track of # consecutively correct (regardless of catch or different)  , but reset
                    % when we get 2 in a row

reversals = 0; % counter for number of changes (reversals) over time
Trial = 0; % counter for each trial 

changeDirection = 1; % Placeholder to remember "direction" of delta changes to count reversals

ordersave   = []; % Saves whether it is order 1 or order 2, which indicates whether the reference frequency is first or second
Modulator   = []; % How much to modulate the dist by
dist_values = []; % Saves the dist values
revValues   = []; % Saves the reversal values, will quit at 12

% Loop until we get the necessary # of reversals or we hit the max number
% of trials allowed

while  reversals < MaxReversals 
    Trial = Trial + 1; % advance the trial every time we go through this loop
    randdiff = 8;
    RandomCatchTrial = randperm(randdiff,1); % 1/8 if the time it is the same and 50% it is different
    
    % Create the comparison stimulus
    if RandomCatchTrial<8
        isDiff =1; 
        n = 2; 
        order = randperm(n, 1); %different trials, order 1 reference is first, order 2 reference is second 
    else
        isDiff=0;
        % In ~1/2 trials the comparison stimulus is same as reference
        order = 3; 
        dist2 = 9999 ; %marks the dist value for the same trials
        dist_values = [dist_values dist2];     
    end;        
    ordersave = [ordersave order];

    for Phase = 1:2 %phase 1 plays first then phase 2
        if (order == 3) || (order == 1 && Phase == 1) || (order == 2 && Phase == 2)
           %if catchtrial, both phases play reference
           %OR if different trial order 1 AND first phase play reference
           %OR if different trial order 2 AND second phase play reference
            
           if order ==1 
               dist_values = [dist_values dist]; %saves all the dist values BEFORE the person hears them/makes a decisions
           end

           soundwav = cos(timeT*2*pi*440); %plays the baseline - doesn't change
        else % makes the comparision stimuli
            if order == 2
                dist_values = [dist_values dist];
            end   %saves all the dist values BEFORE the person hears them/makes a decisions

            newf = 440*2^(dist/12);
            soundwav = cos(timeT*2*pi*newf);
        end

       % Window the beginning and end of the stimuli slightly to avoid clicks
        taper = tukeywin(length(timeT),0.05);
        soundwav = soundwav .* taper';

        if Trial==1 && Phase==1
            disp('Hit return when ready to start') % what the users sees on the screen
            pause
        end
        if Trial~=1 && Phase==1
            disp(['Hit return to begin next trial'])% what the users sees on the screen
            pause
        end

        sound(.05*soundwav*10^(-6/20),fs); %65 dB on MOTU equipment

        if Phase == 1
            pause(2) % Play the two sounds in succession (with a short pause between)
        end

        if Phase == 2
            fprintf('Same (0) or Different (1)?'); % ask them if they were the same or different
            YourAnswer = getkey; %getkey is a function that finds the value of the next keypress, getkey MUST be in the same folder as this code
            fprintf('%c\n', YourAnswer) %displays the key that was entered
            
            if YourAnswer == 49 % key press = 1
                response(Trial) = 1; %the response was 1 (different)
            elseif YourAnswer == 48 % key press = 0
                response(Trial) = 0; %the response was 0 (same)
            end

            % Determine if the user was correct or incorrect
            
            % If it is a catch trial, then isDiff will be 0, and 
            % a response of "Same" (0) will be a match. 
            
            % If it is not a catch trial, isDiff will be 1, and
            % a response of "Different" (1) will be a match.
            
            % In all other cases the trial will not be a match
            match(Trial) = isequal(isDiff,(response(Trial))); 
            
%             outfile2=fullfile( answer{ 1 } ,['Trial_',num2str(Trial)]);
%             save = (outfile2);
          if dist < 0
              dist     = 0.02;
              upstep   = 0.01;
              downstep = upstep / 2.4483; 
          else
              if Trial > 10 %reduce step size
                  upstep = 0.1;
                  downstep = upstep / 2.4483; 
                  if Trial > 20 %reduce step size again (last time)
                      upstep = 0.05;
                      downstep = upstep / 2.4483;
                      if Trial > 30 %reduce step size again (last time)
                          upstep = 0.02;
                          downstep = upstep / 2.4483;
                      end
                  end
              end
          end

        % Adaptively adjust dist as required by 2-down, 1-up
        if (Trial > 1) %nothing happens on trial 1, it does not get included into 2 up and 1 down calculations

            if (match(Trial)) %this is if match = 1, so they were correct this trial

                if (correctInARow == 1) %this means they were correct last trial as well

                    dist = dist - downstep; %they now have got 2 correct in a row, so move the dist down to make it closer to the reference, aka harder

                    if (changeDirection == 1)% if it did change direction count a reversal, record the dist value and then note that it already changed directions
                        reversals = reversals + 1; % counting the reversals
                        revValues = [revValues dist]; %recording the revValues of the dist that JUST happened 
                        changeDirection = -1; %nothing that it just changed directions
                    end;
                    correctInARow = 0; % reset the # in a row counter after you have changed the dist
                else
                    dist = dist; % no change in dist
                    % This is the first correct in a row.  No change in dist
                    correctInARow = 1; % mark our # in a row counter
                end;
            else
                % The subject got the last trial wrong -> increase delta
                dist = dist + upstep;
                if (changeDirection == -1)
                    % This was a reversal in direction of changes
                     reversals = reversals + 1;
                        revValues = [revValues dist];
                        changeDirection = 1;
                end;
                correctInARow = 0;
            end;
        else
            % This was the first trial of the experiments - do nothing
            dist(Trial) = dist(Trial);
        end;

          if dist < 0
                     dist = .02;
                     %disp('it was negative!') 
        end
          if Trial > 84
                    reversals = MaxReversals+1;
          end

    %        outfile2=fullfile( answer{ 1 } ,['Trial_',num2str(Trial)]);
    %         save( outfile2, 'TrialAnswer')

     ResultMatrix = [dist_values', match', response', ordersave'] ;  
     save( outfile , 'ResultMatrix' , 'revValues' ) ;
        end
    
    end
    
  % save( outfile , 'ResultMatrix' , 'revValues', 'ResultMatrixTitle'  ) ; 
end
ans=mean(revValues(end-5:end))
%FullFile =  [ResultMatrixTitle; num2cell(ResultMatrix)];
%save( outfile , 'ResultMatrix' , 'revValues', 'ResultMatrixTitle', 'FullFile'  ) ; 
%ans=mean(revValues(end-6:end));
%FullFile =  [ResultMatrixTitle; num2cell(ResultMatrix)];
%save( outfile , 'ResultMatrix' , 'revValues', 'ResultMatrixTitle', 'FullFile'  ) ;  

% Create a plot of delta against trials
% figure;
% plot(delta,'-o');
% xlabel('Trial number');
% ylabel('\Delta f (Hz)','Interpreter','Tex');

% Calculate the average delta over the last nJNDTrials trials
% JND = mean(delta(end-nJNDTrials:end));
% msgbox(sprintf('JND estimate is %2.2f Hz',JND));