function out = GAP_PitchJND_Passive_500()

close all;
ET = tic;
rng('shuffle');

cue_check = input('Did you check the CueMix Fx configuration is correct? (0 = no, 1 = yes): ');
if cue_check == 0
    errordlg('Check CueMix setting!');
    return
end
 
prompt = {'CueMix MainOut (100 for max, otherwise clock orientation): '}; %This line asks for cuemix trim settings for MainOut. MainOut doesn't have values, so max settings is "100", otherwise, specify clock orientation (for example, for the setting level with the "12" on the bar to the left, it would be 9). 
name = 'MainOut';
numlines = 1;
defaultanswer={'00'}; %input box for cuemix trim setting
answer=inputdlg(prompt, name, numlines, defaultanswer);
cuemixMainout = str2num(answer{1});%main out setting 

prompt = {'Subject ID:',...
          'Session ID:',...
          'Baseline Run:',...
          'Gender ("male" or "female")'};
name = 'Subject Information';
numlines = 1;
defaultanswer = {'null','JNDpitch1','Run3','female'};
answer = inputdlg(prompt, name, numlines, defaultanswer);

if isempty(answer)
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

expParam.subject = answer{1};
expParam.run     = answer{2};
expParam.baseRec = answer{3};
Gender =  answer{4};
project = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
dirs = dfDirs(project);

% create the folder path to save the data files
% baseFilename = ['data\' group, '\', subjectID, '\', session, '\',num_trials,'\'];
dirs.RecFileDir = fullfile(dirs.RecData, expParam.subject, expParam.run);
dirs.SavFileDir = fullfile(dirs.RecData, expParam.subject, expParam.baseRec, [expParam.subject expParam.baseRec 'DRF.mat']);

dirs.tokenDir = fullfile(dirs.RecFileDir, 'speechTokens');
dirs.baseTokenFile = fullfile(dirs.tokenDir,[expParam.subject expParam.run 'BaseToken.wav']);

% check if the foler exists (to avoid overwrite)
if ~exist(dirs.RecFileDir, 'dir')
    mkdir(dirs.RecFileDir);
    while exist(dirs.RecFileDir,'dir') ~= 7
        mkdir(dirs.RecFileDir);
    end
else
    overwrite = inputdlg({'File already exists! Do you want to overwrite?'},'Overwrite',1,{'no'});
    if ~strcmp(overwrite,'yes') & ~strcmp(overwrite,'YES') & ~strcmp(overwrite,'Yes')
        return;
    end
end

if ~exist(dirs.SavFileDir, 'file')
    disp('ERROR: No baseline file at this location!')
    return
end

if ~exist(dirs.tokenDir, 'dir')
    mkdir(dirs.tokenDir);
end

%% Setting up the up-down paradigm (modified based on Palam)
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
UD.cuemixMainout = cuemixMainout;

%% recording audio samples
[BaseToken, fs]= extractSpeechToken(dirs);
subjf0 = calcf0(BaseToken, fs);
PertFreqs = targetf0calc(subjf0, UD.xMax, UD.xMin);
numPertFreqs = length(PertFreqs);
PertTokens = generatef0JNDTokens(dirs, numPertFreqs, PertFreqs);

%%%%%Visual Presentation
[h2, h3, h4] = JNDVisualPresentation;
pause(5);

tr = 0;
while (UD.stop == 0) & tr < UD.totalTrials
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %HERE IS THE MAGIC!!!!
    sound(Token1, fs)
    pause(TokenLen1 + 0.01 + UD.ISI)
    sound(Token2, fs)
    pause(TokenLen2 + 0.01)
    %HERE IS ALL YOU HAVE BEEN WAITING FOR!!! 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Present the YES/NO question
    %     set(h1,'Visible','off');
    set(h2, 'String','PITCH','FontSize',80)%was 'WERE THEY DIFFERENT'
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
    set(h2, 'String','','FontSize',120)
    set(h3, 'Visible','off');
    set(h4, 'Visible','off');
    drawnow
    if conVar == 1 % update the UD structure on real trials
        UD = adaptiveUD_Update(UD, response);
        UD.catchResponse(tr,1) = NaN;
    else % when it is a catch trial do not update UD structure (i.e., do not change the up-down steps based on catch trials)
        UD.catchResponse(tr,1) = response;
    end
    pause(1) %this is between two trials   
end
close all;
elapsed_time = toc(ET)/60;
disp (sprintf('Total time: %f (min)',elapsed_time));

UD.BaseToken = BaseToken;
UD.PertTokens = PertTokens;
%CES - To disable saving of reactiontimes (we are doing this to avoid confusion since the investigator is keying in the selection at present)
UD.reactionTime = ones(size(ReactionTime))*10000;

expFiles = fullfile(dirs.RecFileDir, 'ExperimentalParameters.mat');
save(expFiles, 'UD');

dataFileName = fullfile(dirs.RecFileDir, [expParam.subject expParam.run 'data.mat']);
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
        save(dataFileName, 'dataFile');
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

function f0 = calcf0(x,fs)
% Created by Gabriel Galindo
% Formatted by Dante Smith -12/11/15

lim_inf = ceil(fs/(500));
lim_sup = floor(fs/(50));
U = xcov(x,'unbias');
U = U(ceil(end/2):end);
U = (U(lim_inf:lim_sup)-min(U(lim_inf:lim_sup)))/(max(U(lim_inf:lim_sup)) - min(U(lim_inf:lim_sup)));
[M,P] = findpeaks(U);

if isempty(P)
    f0 = NaN;
else
    P = P(find(M >= 0.9,1,'first'));
    if isempty(P)
        f0 = NaN;
    else
        f0 = fs/(P + lim_inf);
    end

    NFFT = pow2(nextpow2(length(x)/4));
    [Pxx,Fxx] = pwelch(x,NFFT,[],[],fs,'onesided');

    if ~isnan(f0)
        H = Pxx(find(Fxx>=f0,1,'first'));
        if (10*log10(max(Pxx)/H) > 80)
            f0 = NaN;
        end
    end   
end
end

function freqs = targetf0calc(f0, maxC, minC)
%calculates all possible freq vals spaced 1 cent apart. 

numCents = maxC - minC;

for i = 1:numCents
    freqs(i) = f0*2^(i/1200);
end
end