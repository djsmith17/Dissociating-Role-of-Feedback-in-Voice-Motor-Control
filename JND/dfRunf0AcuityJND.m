function dfRunf0AcuityJND()
% Received from Ayoub Daliri
close all;
ET = tic;
rng('shuffle');

prompt = {'Subject ID:',...
          'Session ID:',...
          'Baseline Run:',...
          'Baseline Trial:',...
          'Direction ("Above" or "Below"):',...
          'Gender ("male" or "female")'};
name = 'Subject Information';
numlines = 1;
defaultanswer = {'null','fA1', 'BV1', '3', 'Above', 'female'};
answer = inputdlg(prompt, name, numlines, defaultanswer);

if isempty(answer)
    return
end

num_trials = questdlg('Practice or Full?','Length','Practice','Full','Full') ;
switch num_trials
    case 'Practice'
        totalTrials = 9;
    case 'Full'
        totalTrials = 60; %max number of trials if max trials/reversals not reached
end

UD.project = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
UD.subject = answer{1};
UD.run     = answer{2};
UD.baseRec = answer{3};
UD.baseTrial = str2double(answer{4});
UD.JNDDirection = answer{5};
UD.gender =  answer{6};

dirs = dfDirs(UD.project);
% Folder paths to save data files
dirs.RecFileDir = fullfile(dirs.RecData, UD.subject, UD.run);
dirs.SavFileDir = fullfile(dirs.RecData, UD.subject, UD.baseRec, [UD.subject UD.baseRec 'DRF.mat']);

dirs.tokenDir = fullfile(dirs.RecFileDir, 'speechTokens');
dirs.baseTokenFile = fullfile(dirs.tokenDir,[UD.subject UD.run 'BaseToken.wav']);

if ~exist(dirs.RecFileDir, 'dir')
    mkdir(dirs.RecFileDir);
end

if ~exist(dirs.SavFileDir, 'file')
    disp('ERROR: No voice file at this location!')
    return
end

if exist(dirs.tokenDir, 'dir')
    rmdir(dirs.tokenDir, 's')  
end
mkdir(dirs.tokenDir);

% Setting up the up-down paradigm (modified based on Palam)
UD.totalTrials = totalTrials;
UD.up = 1;    % Number of consecutive responses before an increase
UD.down = 2;  % Number of consecutive responses before a decrease
stepSize = 4; %This is something to tune; in cents
UD.stepSizeUp = stepSize; %Levitt (1971) 2/1 rule for 71% in MacMillian Chapter 11 with step per Garcia-Perez (1998); Was: Size of step up ; stepSize/ .5488 ensures 80.35 % correct; see Garcia-Perez 1998
UD.stepSizeDown = stepSize; % Size of step down
UD.stopCriterion = 'reversals'; % stop the procedure based on number of 'trials' | 'reversals'
UD.stopRule = 10;  %stop procedure after this number of trials/reversals
UD.startValue = 50; % initial difference in cents between speaker's fo and fo of stimulus in headphones
UD.xMax = 200; %max difference between speaker's fo and fo of stimulus in headphones
UD.xMin = 0; %min difference between speaker's fo and fo of stimulus in headphones
UD.xAll = UD.xMin:UD.xMax;
UD.xLen = length(UD.xAll);
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
waitForKeyPress = 3 ; % in seconds
UD.ISI = .5; %Interstimulus interval (ISI) within each trial (between stimulus 1 and stimulus 2 in pair) in seconds
UD.measuredDelay = 0.0; %Measured delay of instruments to be incoportated for accurate ISI and token length

if strcmp(UD.JNDDirection, 'Above')
    sign = 1;
else
    sign = -1;
end

% Generate audio tokens
[BaseToken, fs]= dfGenerateBT(dirs, UD.baseTrial); %Extract a Speech Token. Located in JND Folder
subjf0 = dfcalcf0Praat(dirs);                      %Calculate f0 using praat. Located in JND Folder
PertFreqs = targetf0calc(subjf0, UD.xMax, UD.xMin, sign); %Located Below
numPertFreqs = length(PertFreqs);
PertTokens = dfGeneratePT(dirs, numPertFreqs, PertFreqs); %Generate Pert Tokens. Located in JND Folder

fprintf('Starting f0 Acuity Task for %s with f0 of %f\n', UD.subject, subjf0)
%%%%%Visual Presentation
[h2, h3, h4] = JNDVisualPresentation;
pause(5);

tr = 0;
while (UD.stop == 0) && tr < UD.totalTrials
    tr = tr +1;
    %Present the word
    set(h2,'String','+')
    drawnow;
    
    tempVar = randperm(5);
    Pert = UD.xCurrent; %cents
    PertToken = PertTokens(Pert, :);
    if tempVar(1) == 1 || tempVar(1) == 3   % scenario I (first one is Pert) : % 40% of trials
        Token1 = PertToken;
        Token2 = BaseToken;
        conVar = 1;
    elseif tempVar(1) == 2 || tempVar(1) == 4 % scenario II (first one is no Pert) : % 40% of trials
        Token1 = BaseToken;
        Token2 = PertToken;
        conVar = 1;
    else % Catch trials : 20% of trials
        Token1 = BaseToken;      
        Token2 = BaseToken;
        conVar = 0; %catch trials will be randomly presented but will not be included in the adaptive procedure        
    end
    TokenLen1 = length(Token1)/fs; TokenLen2 = length(Token2)/fs;
    
    JNDMessage(tr, Pert, conVar, 0, 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %HERE IS THE MAGIC!!!!
    sound(Token1, fs)
    pause(TokenLen1 + UD.ISI + UD.measuredDelay)
    sound(Token2, fs)
    pause(TokenLen2 + UD.measuredDelay)
    %HERE IS ALL YOU HAVE BEEN WAITING FOR!!! 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Present the YES/NO question
    set(h2, 'String', 'PITCH', 'FontSize', 80)
    set(h3, 'Visible','on');
    set(h4, 'Visible','on');
    drawnow
    keyCorrect = 1;
    while keyCorrect
        
        keyCorrect = 0;
        [bb, ReactionTime(tr)] = GetKey_Ayoub(1, waitForKeyPress);
        
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
    JNDMessage(tr, Pert, conVar, response, 2);
    
    set(h2, 'String','','FontSize',120)
    set(h3, 'Visible','off');
    set(h4, 'Visible','off');
    drawnow
    if conVar == 1 % update the UD structure on real trials
        UD = dfAdaptiveUpdateJND(UD, response);
        UD.catchResponse(tr,1) = NaN;
    else % when it is a catch trial do not update UD structure (i.e., do not change the up-down steps based on catch trials)
        UD.catchResponse(tr,1) = response;
    end
    
    
    pause(1) %this is between two trials   
end
close all;
elapsed_time = toc(ET)/60;
fprintf('\nElapsed Time: %f (min)\n', elapsed_time)

UD.subjf0    = subjf0;
UD.BaseToken = BaseToken;
UD.PertTokens = PertTokens;
UD.reactionTime = ones(size(ReactionTime))*10000;
UD.elapsedTime = elapsed_time;

UD.performedTrials = length(UD.catchResponse);
UD.JNDTrials = length(UD.reversal);
UD.catchTrials = length(UD.catchResponse) - length(UD.reversal);
UD.reversals = max(UD.reversal);
UD.catchCorrect = sum(UD.catchResponse == 0);
UD.catchAccuracy = 100*(UD.catchCorrect/UD.catchTrials);

expFiles = fullfile(dirs.RecFileDir, [UD.subject UD.run 'DRF.mat']);

switch num_trials
    case 'Practice'
        return
    case 'Full'
        save(expFiles, 'UD'); %Only save if it was a full set of trials
end

end

function [h2, h3, h4] = JNDVisualPresentation
monitorSize = get(0,'Monitor');
if size(monitorSize,1) == 1
    figPosition = [1 200 monitorSize(3) monitorSize(4)-200];
elseif size(monitorSize,1) == 2
    figPosition = [monitorSize(2,1) monitorSize(2,2) monitorSize(1,3) monitorSize(2,4)];
end

figure1 = figure('Color',[0 0 0],'Menubar','none','Position', figPosition);

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
end

function freqs = targetf0calc(f0, maxC, minC, sign)
%calculates all possible freq vals spaced 1 cent apart. 

numCents = maxC - minC;

for i = 1:numCents
    freqs(i) = f0*2^(sign*i/1200);
end
end

function JNDMessage(tr, Pert, conVar, response, state)

if state == 1
    msg = ['Trial ' num2str(tr) ' at ' num2str(Pert) ' cents: '];
    
    if conVar == 1
        msg = [msg 'Is Diff, Answered '];
    else
        msg = [msg 'Is Same, Answered '];
    end
else
    if response == 1
        msg = 'Diff\n';
    else
        msg = 'Same\n';
    end
end

fprintf(msg)
end