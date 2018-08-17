function dfRunf0AcuityJND_AXB()
%dfRunf0Acuity_AXB() is an experimental paradigm for testing perceptual
%acuity to pitch using a AXB method. This program uses a previously
%generated set of speech tokens to be presented to the participant. On each
%trial of this task, the participant listens to three tokens of their own
%voice, and must decide if the pitch of the middle token is different than
%the pitch of the first token, or the pitch of the last token. Left and
%Right arrow keys are used as input, and there is visual presentation of
%the task. The task adaptively changes to find a threshold at which
%participants can discriminate between two pitches. Hitting the Esc key 
%after a token presentation will exit the script. The core version of the
%script was received from Ayoub Daliri.
%
% The Voice Token Data file will load a structure that is named GT
%
% This script is dependent on the following external functions:
% -dfDirs
% -dfAdaptiveUpdateJNDAXB.m
% -GetKey_Ayoub.m
%
% This script has the following subfunctions:
% -JNDVisualPresentation
% -JNDMessage
% -pseudoRandomTrialOrder
% -accuLogic

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
defaultanswer = {'null', 'fAX1', 'GT1', 'Diff', 'female'};
answer = inputdlg(prompt, name, numlines, defaultanswer);

if isempty(answer)
    return
end

num_trials = questdlg('Practice or Full?','Length','Practice','Full','Full') ;
switch num_trials
    case 'Practice'
        totalTrials = 20;
    case 'Full'
        totalTrials = 100; % Max number of trials if # reversals not reached
end

UD.project   = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
UD.subject   = answer{1};
UD.run       = answer{2};
UD.tokenFile = answer{3};
UD.inst      = answer{4};
UD.gender    = answer{5};

dirs = dfDirs(UD.project);
% Folder paths to save data files
dirs.RecFileDir = fullfile(dirs.RecData, UD.subject, UD.run);
dirs.TokenFile  = fullfile(dirs.RecData, UD.subject, UD.tokenFile, [UD.subject UD.tokenFile 'DRF.mat']);

if ~exist(dirs.RecFileDir, 'dir')
    mkdir(dirs.RecFileDir);
end

if ~exist(dirs.TokenFile, 'file')
    disp('ERROR: No tokens at this location! Please try another')
    return
else
    load(dirs.TokenFile); % Returns a struct 'GT'.
end

%Token Generation Output;
UD.baseRec    = GT.baseRec;
UD.baseTrial  = GT.baseTrial;
UD.subjf0     = GT.subjf0;
UD.PertFreqs  = GT.PertFreqs;
UD.fs         = GT.fs;
UD.BaseToken  = GT.BaseToken;
UD.PertTokens = GT.PertTokens;

% Setting up the up-down paradigm (modified based on Palam)
UD.totalTrials = totalTrials;
UD.up = 1;    % Number of consecutive responses before an increase
UD.down = 2;  % Number of consecutive responses before a decrease
stepSize = 4; % This is something to tune; in cents
UD.stepSizeUp = stepSize; %Levitt (1971) 2/1 rule for 71% in MacMillian Chapter 11 with step per Garcia-Perez (1998); Was: Size of step up ; stepSize/ .5488 ensures 80.35 % correct; see Garcia-Perez 1998
UD.stepSizeDown = stepSize; % Size of step down
UD.BIGstep      = 10;
UD.smallStep    = 1;
UD.stopCriterion = 'reversals'; % stop the procedure based on number of 'trials' | 'reversals'
UD.stopRule      = 10; % stop procedure after this number of trials/reversals
UD.startValue    = 50; % initial difference in cents between speaker's fo and fo of stimulus in headphones
UD.xMax = 200;         % max difference between speaker's fo and fo of stimulus in headphones
UD.xMin = 1;           % min difference between speaker's fo and fo of stimulus in headphones
UD.xAll = -100:0.5:100;
UD.xAll = UD.xAll(~logical(UD.xAll == 0));
UD.xLen = length(UD.xAll);
UD.truncate = 'yes';
UD.response = [];
UD.stop = 0;
UD.u = 0;
UD.d = 0;
UD.direction = [];
UD.reversal  = 0;
UD.xCurrent  = UD.startValue;
UD.x = UD.startValue;
UD.xStaircase = [];
UD.catchResponse = [];
UD.allTrialPerts = [];
UD.allTrialTypes = [];
waitForKeyPress  = 3;   % Seconds
UD.ISI           = 0.5; % Interstimulus interval (ISI) within each trial (between stimulus 1 and stimulus 2 in pair) in seconds
UD.measuredDelay = 0.0; % Measured delay of instruments to be incoportated for accurate ISI and token length

fprintf('Starting f0 Acuity Task for %s with f0 of %f\n\n', UD.subject, UD.subjf0)
%%%%%Visual Presentation
[ctrMsg, leftArw, righArw] = JNDVisualPresentation;
pause(5);

tr = 0;
countSD = [0 0];
while (UD.stop == 0) && tr < UD.totalTrials
    tr = tr + 1;
    % Present the word
    set(ctrMsg, 'String', '+')
    drawnow;
       
    PertDist  = UD.xCurrent; %cents
    hPertDist = PertDist/2;
    
    [trialInd, countSD] = pseudoRandomTrialOrder(4, countSD);
    if trialInd == 1       % A and X Match, B is different
        pertA  = hPertDist;
        pertX  = hPertDist;
        pertB  = -hPertDist;        
        conVar = 1;
    elseif trialInd == 3   % A and X Match, B is different
        pertA  = -hPertDist;
        pertX  = -hPertDist;
        pertB  = hPertDist;
        conVar = 1;
    elseif trialInd == 2   % B and X Match, A is different
        pertA  = hPertDist;
        pertX  = -hPertDist;
        pertB  = -hPertDist;        
        conVar = 0;
    elseif trialInd == 4   % B and X Match, A is different
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
    set(ctrMsg, 'String', 'PITCH', 'FontSize', 80)
    set([leftArw righArw], 'Visible','on');
    drawnow
    keyCorrect = 1;
    while keyCorrect
        
        keyCorrect = 0;
        [bb, ReactionTime(tr)] = GetKey_Ayoub(1, waitForKeyPress);
        
        % wait until a correct key is pressed
        if (bb ~= 28) & (bb ~= 29) & (bb ~= 27)
            keyCorrect = 1;
        end
        if isempty(bb)
            keyCorrect = 1;
        end
        
        if bb == 28            %28 is "FIRST" | Left ARROW KEY; 
            response = 1;
        elseif bb == 29        %29 is "LAST"  | Right ARROW KEY;
            response = 0;
        elseif bb == 27        %27 is ESC     | Quit for now and save;
            response = 0;
            UD.stop = 1; 
        end      
    end
    JNDMessage(tr, trialPerts, PertDist, conVar, response, 2);  
    
    set(ctrMsg, 'String', '', 'FontSize',120)
    set([leftArw righArw], 'Visible','off');
    drawnow 
    
    [trialType, correct] = accuLogic(UD.inst, conVar, response);
     
    UD = dfAdaptiveUpdateJNDAXB(UD, response, correct);
    UD.catchResponse(tr,1) = correct;
    UD.allTrialPerts = cat(1, UD.allTrialPerts, trialPerts);
    UD.allTrialTypes = cat(1, UD.allTrialTypes, trialType);
    pause(1) %Inter-Trial time  
end
close all;
elapsed_time = toc(ET)/60;
fprintf('\nElapsed Time: %f (min)\n', elapsed_time)

UD.reactionTime = ones(size(ReactionTime))*10000;
UD.elapsedTime  = elapsed_time;

UD.performedTrials = length(UD.catchResponse);
UD.JNDTrials       = length(UD.reversal);
UD.catchTrials     = sum(~isnan(UD.catchResponse));
UD.reversals       = max(UD.reversal);
UD.catchCorrect    = sum(UD.catchResponse);
UD.catchAccuracy   = 100*(UD.catchCorrect/UD.catchTrials);

expFiles = fullfile(dirs.RecFileDir, [UD.subject UD.run 'DRF.mat']);
switch num_trials
    case 'Practice'
        return
    case 'Full'
        fprintf('Saving JND recording\n')
        save(expFiles, 'UD'); % Only save if it was a full set of trials
end
end

function [ctrMsg, leftArw, righArw] = JNDVisualPresentation
% JNDVisualPresentation generates a figure on the computer monitor to 
% provide trial information. Some of these positions are a little hard 
% baked. In future versions, I will make them modular to the screen size.
%
% If there are two monitors, this will scale the window to the size of the
% second screen. If there is only one monitor, it will create a small
% window, with the thought that it makes things easier for testing on a
% single-screen computer.

monitorSize = get(0,'Monitor');
numMon = size(monitorSize, 1);

if numMon == 1
    W = monitorSize(3);
    H = monitorSize(4);
    W2 = W/2;
    H2 = H/2;
    XPos = W2;
    YPos = 50;
    
    figPosition = [XPos YPos W2 H2];
elseif numMon == 2
    [~, mon] = max(monitorSize(:,1));
    
    figPosition = [monitorSize(mon,1) monitorSize(mon,2) monitorSize(1,3) monitorSize(1,4)];
end
figDim.winPos = figPosition;

% Ready Annotation Dim
figDim.rdAnoD = [700 300];
figDim.rdAnoW = round(figDim.rdAnoD(1)/figDim.winPos(3), 2); 
figDim.rdAnoH = round(figDim.rdAnoD(2)/figDim.winPos(4), 2);
figDim.rdAnoX = 0.5 - figDim.rdAnoW/2;
figDim.rdAnoY = 0.5 - figDim.rdAnoH/2;
figDim.rdAnoPos = [figDim.rdAnoX figDim.rdAnoY figDim.rdAnoW figDim.rdAnoH];

figure1 = figure('Color', [0 0 0], 'Position', figDim.winPos, 'MenuBar', 'none');

ctrMsg = annotation(figure1, 'textbox', figDim.rdAnoPos,...
                             'Color', [1 1 1],...
                             'String', 'READY',...
                             'LineStyle', 'none',...
                             'HorizontalAlignment', 'center',...
                             'VerticalAlignment', 'middle',...
                             'FontSize', 130,...
                             'FontName', 'Arial',...
                             'FitBoxToText', 'off',...
                             'EdgeColor', 'none',...
                             'BackgroundColor', [0 0 0],...
                             'Visible','on');

leftArw = annotation(figure1, 'textbox', [0.025 0.15 0.45 0.3],...
                              'Color', [0 0 0],...
                              'String', {'< FIRST'},...
                              'HorizontalAlignment', 'center',...
                              'VerticalAlignment', 'middle',...
                              'FontSize', 60,...
                              'FontName', 'Arial',...
                              'LineStyle', 'none',...
                              'BackgroundColor', [1 1 1],...
                              'Visible','off');

righArw = annotation(figure1, 'textbox', [0.52 0.15 0.45 0.3],...
                              'Color',[0 0 0],...
                              'String', {'LAST >'},... 
                              'HorizontalAlignment', 'center',...
                              'VerticalAlignment', 'middle',...
                              'FontSize', 60,...
                              'FontName', 'Arial',...
                              'LineStyle', 'none',...
                              'BackgroundColor', [1 1 1],...
                              'Visible','off');

drawnow;
end

function JNDMessage(tr, PertVals, PertDist, conVar, response, state)
% JNDMessage() provides output to the researcher about the current trial
% taking place, and the participant's live response. It is useful to have 
% this information while testing JND to confirm that the logic of the
% staircase adaptive JND is working correctly. Accuracy, and reversals are
% not currently tracked here, but can be derived implictely by observation.

if state == 1 % Token Presentation
    msg = ['Trial ' num2str(tr) ' at ' num2str(PertDist) ' (' num2str(PertVals(1)) ', ' num2str(PertVals(2)) ', ' num2str(PertVals(3)) '): '];
    
    if conVar == 1 % A and X Match, B is different
        msg = [msg 'Is Last, Answered '];
    else           % B and X Match, A is different
        msg = [msg 'Is First, Answered '];
    end
else          % Subject's Response
    if response == 1
        msg = 'First\n'; % Answered First
    else
        msg = 'Last\n';  % Answered Last
    end
end

fprintf(msg)
end

function [trialInd, countSD] = pseudoRandomTrialOrder(nTrialTypes, countSD)
% pseudoRandomTrialOrder(nTrialTypes, countSD) ensures that there are no
% more than the max number of trials of a type in row (e.g. AAA or BBB)

maxCount = 3; % Currently hard set to 3 in a row

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

function [trialType, correct] = accuLogic(inst, conVar, response)
% accuLogic(UD, conVar, response) determines if the response by the
% participant was correct or incorrect. It has outputs for two instructions
% types: Whether the subject was asked to respond for the token that was 
% the 'same' as 'X', or if they were asked to respond for the token that 
% was 'different' than 'X'. If we finalize an instruction, this function
% can be made much simplier.
%
% INPUTS:
% inst:     Instruction. Either 'Same' or 'Diff'
% conVar:   Condition. 1 = A and X Match, B is different. 
%                      0 = B and X Match, A is different.
% response: Subject response. 1 = 'First'.
%                             0 = 'Last'.
%
% OUTPUTS:
% trialType: Result of trial condition and correctness. This is used in
%            conjunction with the other output variable 'correct' to better
%            describe which trials were correct/incorrect to observe any 
%            bias. 
% correct:   correct (1) or incorrect (0)

if strcmp(inst, 'Same')
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
