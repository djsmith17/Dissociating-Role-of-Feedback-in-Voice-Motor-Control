function dfRunf0AcuityJND_AXB()
%dfRunf0Acuity_AXB() is an experimental paradigm for testing perceptual
%acuity to pitch using a AXB method. This program uses a previously
%generated set of speech tokens to be presented to the participant. On each
%trial of this task, the participant listens to three tokens of their own
%voice, and must decide if the pitch of the middle token is different than
%the pitch of the first token, or the pitch of the last token. Left and
%Right arrow keys are used as input, and there is visual presentation of
%the task. The task adaptively changes to find a threshold at which
%participants can discriminate between two pitches. The core version of the
%script was received from Ayoub Daliri.
%
%The Voice Token Data file will load a structure that is called GT
%
%This script makes use of the following outside functions:
%-dfAdaptiveUpdateJNDAXB.m
%-GetKey_Ayoub.m

close all;
ET = tic;
rng('shuffle');

prompt = {'Subject ID:',...
          'Session ID:',...
          'Tokens File:',...
          'Instruction ("Same" or "Diff):',...
          'Gender ("male" or "female")'};
name = 'Subject Information';
numlines = 1;
defaultanswer = {'null','fAX1', 'GT1', 'Diff', 'female'};
answer = inputdlg(prompt, name, numlines, defaultanswer);

if isempty(answer)
    return
end

num_trials = questdlg('Practice or Full?','Length','Practice','Full','Full') ;
switch num_trials
    case 'Practice'
        totalTrials = 20;
    case 'Full'
        totalTrials = 100; %max number of trials if max trials/reversals not reached
end

UD.project = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
UD.subject = answer{1};
UD.run     = answer{2};
UD.tokenFile = answer{3};
UD.inst    = answer{4};
UD.gender  = answer{5};

dirs = dfDirs(UD.project);
% Folder paths to save data files
dirs.RecFileDir = fullfile(dirs.RecData, UD.subject, UD.run);
dirs.TokenFile  = fullfile(dirs.RecData, UD.subject, UD.tokenFile, [UD.subject UD.tokenFile 'DRF.mat']);

if ~exist(dirs.RecFileDir, 'dir')
    mkdir(dirs.RecFileDir);
end

if ~exist(dirs.TokenFile, 'file')
    disp('ERROR: No tokens at this location!')
    return
end

%Token Generation Output;
load(dirs.TokenFile); %Should return struct GT. This should become a try statement
UD.baseRec    = GT.baseRec;
UD.baseTrial  = GT.baseTrial;
UD.subjf0     = GT.subjf0;
UD.pertFreqs  = GT.pertFreqs;
UD.fs         = GT.fs;
UD.BaseToken  = GT.BaseToken;
UD.PertTokens = GT.PertTokens;

% Setting up the up-down paradigm (modified based on Palam)
UD.totalTrials = totalTrials;
UD.up = 1;    % Number of consecutive responses before an increase
UD.down = 2;  % Number of consecutive responses before a decrease
stepSize = 4; %This is something to tune; in cents
UD.stepSizeUp = stepSize; %Levitt (1971) 2/1 rule for 71% in MacMillian Chapter 11 with step per Garcia-Perez (1998); Was: Size of step up ; stepSize/ .5488 ensures 80.35 % correct; see Garcia-Perez 1998
UD.stepSizeDown = stepSize; % Size of step down
UD.BIGstep      = 10;
UD.smallStep    = 1;
UD.stopCriterion = 'reversals'; % stop the procedure based on number of 'trials' | 'reversals'
UD.stopRule      = 10; % stop procedure after this number of trials/reversals
UD.startValue    = 50; % initial difference in cents between speaker's fo and fo of stimulus in headphones
UD.xMax = 200; %max difference between speaker's fo and fo of stimulus in headphones
UD.xMin = 1; %min difference between speaker's fo and fo of stimulus in headphones
UD.xAll = -100:0.5:100;
UD.xAll = UD.xAll(~logical(UD.xAll == 0));
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
UD.catchResponse = [];
UD.allTrialPerts = [];
UD.allTrialTypes = [];
waitForKeyPress = 3 ; % in seconds
UD.ISI = .5; %Interstimulus interval (ISI) within each trial (between stimulus 1 and stimulus 2 in pair) in seconds
UD.measuredDelay = 0.0; %Measured delay of instruments to be incoportated for accurate ISI and token length

fprintf('Starting f0 Acuity Task for %s with f0 of %f\n\n', UD.subject, subjf0)
%%%%%Visual Presentation
[h2, h3, h4] = JNDVisualPresentation;
pause(5);

tr = 0;
countSD = [0 0];
while (UD.stop == 0) && tr < UD.totalTrials
    tr = tr + 1;
    %Present the word
    set(h2,'String','+')
    drawnow;
       
    PertDist = UD.xCurrent; %cents
    hPertDist = PertDist/2;
    
    [trialInd, countSD] = pseudoRandomTrialOrder(4, countSD);
    if trialInd == 1       %Matches A
        pertA  = hPertDist;
        pertX  = hPertDist;
        pertB  = -hPertDist;        
        conVar = 1;
    elseif trialInd == 3   %Matches A
        pertA  = -hPertDist;
        pertX  = -hPertDist;
        pertB  = hPertDist;
        conVar = 1;
    elseif trialInd == 2   %Matches B
        pertA  = hPertDist;
        pertX  = -hPertDist;
        pertB  = -hPertDist;        
        conVar = 0;
    elseif trialInd == 4   %Matches B
        pertA  = -hPertDist;
        pertX  = hPertDist;
        pertB  = hPertDist;
        conVar = 0;
    end
    trialPerts = [pertA pertX pertB];
    
    indA = find(UD.xAll == pertA); indX = find(UD.xAll == pertX); indB = find(UD.xAll == pertB);
    TokenA = UD.PertTokens(indA, :);  TokenX = UD.PertTokens(indX, :);  TokenB = UD.PertTokens(indB, :);
    TokenLenA = length(TokenA)/UD.fs; TokenLenX = length(TokenX)/UD.fs; TokenLenB = length(TokenB)/UD.fs;
    
    JNDMessage(tr, trialPerts, PertDist, conVar, 0, 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %HERE IS THE MAGIC!!!!
    sound(TokenA, UD.fs)
    pause(TokenLenA + UD.ISI + UD.measuredDelay)
    sound(TokenX, UD.fs)
    pause(TokenLenX + UD.ISI + UD.measuredDelay)
    sound(TokenB, UD.fs)
    pause(TokenLenB + UD.measuredDelay)
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
        
        if bb == 28            %28 is "FIRST" | Left ARROW KEY; 
            response = 1;
        elseif bb == 29        %29 is "LAST"  | Right ARROW KEY;
            response = 0;
        end      
    end
    JNDMessage(tr, trialPerts, PertDist, conVar, response, 2);  
    
    set(h2, 'String','','FontSize',120)
    set(h3, 'Visible','off');
    set(h4, 'Visible','off');
    drawnow 
    
    [trialType, correct] = accuLogic(UD, conVar, response);
     
    UD = dfAdaptiveUpdateJNDAXB(UD, response, correct);
    UD.catchResponse(tr,1) = correct;
    UD.allTrialPerts = cat(1, UD.allTrialPerts, trialPerts);
    UD.allTrialTypes = cat(1, UD.allTrialTypes, trialType);
    pause(1) %this is between two trials   
end
close all;
elapsed_time = toc(ET)/60;
fprintf('\nElapsed Time: %f (min)\n', elapsed_time)

UD.reactionTime = ones(size(ReactionTime))*10000;
UD.elapsedTime  = elapsed_time;

UD.performedTrials = length(UD.catchResponse);
UD.JNDTrials = length(UD.reversal);
UD.catchTrials = sum(~isnan(UD.catchResponse));
UD.reversals = max(UD.reversal);
UD.catchCorrect = sum(UD.catchResponse);
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
    'String',{'< FIRST'},...
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
    'String',{'LAST >'},... 
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

function JNDMessage(tr, PertVals, PertDist, conVar, response, state)

if state == 1
    msg = ['Trial ' num2str(tr) ' at ' num2str(PertDist) ' (' num2str(PertVals(1)) ', ' num2str(PertVals(2)) ', ' num2str(PertVals(3)) '): '];
    
    if conVar == 1
        msg = [msg 'Is Last, Answered '];
    else
        msg = [msg 'Is First, Answered '];
    end
else
    if response == 1
        msg = 'First\n';
    else
        msg = 'Last\n';
    end
end

fprintf(msg)
end

function [trialInd, countSD] = pseudoRandomTrialOrder(nTrialTypes, countSD)
maxCount = 3;

tempVar = randperm(nTrialTypes);
trialInd = tempVar(1);
if trialInd == 1 || trialInd == 3     % A Trials
    countSD(1) = 0;
    countSD(2) = countSD(2) + 1;
    if countSD(2) == (maxCount + 1)
        countSD(1) = countSD(1) + 1;
        countSD(2) = 0;
        trialInd = 2*round(rand) + 2;
    end           
elseif trialInd == 2 || trialInd == 4 % B Trials
    countSD(1) = countSD(1) + 1;
    countSD(2) = 0;
    if countSD(1) == (maxCount + 1)
        countSD(1) = 0;
        countSD(2) = countSD(2) + 1;
        trialInd = 2*round(rand) + 1;
    end     
end
end

function [trialType, correct] = accuLogic(UD, conVar, response)

if strcmp(UD.inst, 'Same')
    if conVar == 1
        if response == 1 
            trialType = 1; %Correct AX
            correct  = 1;  %Correct
        else
            trialType = 2; %Incorrect AX
            correct  = 0;  %Incorrect
        end
    elseif conVar == 0        
        if response == 0 
            trialType = 3; %Correct XB
            correct  = 1;  %Correct
        else
            trialType = 4; %Incorrect XB
            correct  = 0;  %Incorrect
        end
    end
else
    if conVar == 1
        if response == 0 
            trialType = 1; %Correct AX
            correct  = 1;  %Correct
        else
            trialType = 2; %Incorrect AX
            correct  = 0;  %Incorrect
        end
    elseif conVar == 0        
        if response == 1 
            trialType = 3; %Correct XB
            correct  = 1;  %Correct
        else
            trialType = 4; %Incorrect XB
            correct  = 0;  %Incorrect
        end
    end
end 
end
