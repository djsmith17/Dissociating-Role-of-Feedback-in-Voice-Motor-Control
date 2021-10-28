function VisFeedback_testbed
% VisFeedback_testbed() is a script which creates an environment to test 
% the elements created by 'presentation\dfSetVisFB()' and updated
% using 'presentation\dfUpdateVisFB()'.
%
% The window is a bare-bones example of how the visualizations and
% instructions will be presented to the participant, without any additional
% data recording, or behavioral testing requirements. 
% This script includes sub-functions which outline the trial progression
% and organize the presentation/removal of visualization elements.
%
% One feature that is especially highlighted in this testbed is 
% the ability to provide feedback to the participant about their 
% loudness or vocal intensity. For more information about how to 
% adjust these elements and others, please refer to the critical 
% functions:
% -presentation\dfSetVisFB()
% -presentation\dfUpdateVisFB()
%
% Author: Dante J Smith
% Updated: 10/28/2021

% Current Recording Session
curSess     = 'VisFBTest';

% Which monitor to use (1 or 2)
defMon = 2;

% Random vocalization vector
recordedRMS = (5+rand(1991,1))/1000;
meanRMS     = mean(recordedRMS);

% Num of trials to perform, if relevant
numTrial = 4;

% Vocalization intensity parameters
targRMS   = 50; % Participant vocalization intensity target (dB)
boundsRMS = 3;  % Acceptable upper and lower bounds for target itensity (+/- dB)

close all
[msrStr, annoStr] = dfSetVisFB(defMon, curSess, targRMS, boundsRMS);

% Paradigm examples to test 
pertParadigm(msrStr, annoStr, meanRMS)  % Testing perturbation paradigm
% triggerTest(annoStr, numTrial)         % Testing Optical Trigger
% endoParadigm(annoStr, numTrial)        % Testing Endoscopy paradigm

close all
end

function pertParadigm(msrStr, annoStr, meanRMS)
% This sub-function mimcs the ordering of visual feedback used
% for experiments studying unexpected perturbations to sensory 
% feedback

% Close the curtains
pause(5); % Let them breathe a sec
set(annoStr.Ready,'Visible','off');

% 
set(annoStr.plus,'Visible','on');
pause(1.0)
set(annoStr.plus,'Visible','off');

set([annoStr.EEE annoStr.visTrig], 'Visible','on');
pause(4.0)
set([annoStr.EEE annoStr.visTrig], 'Visible','off');

[color, newPos, loudResult] = dfUpdateVisFB(msrStr, meanRMS);

set(annoStr.LoudRec, 'position', newPos);
set(annoStr.LoudRec, 'Color', color, 'FaceColor', color);

set([annoStr.LoudRec, annoStr.fbLines], 'Visible','on'); 
pause(2.0)
set([annoStr.LoudRec, annoStr.fbLines], 'Visible','off');
end

function triggerTest(annoStr, numTrial)
% This sub-function mimcs the presentation of a optical trigger
% in the lower left corner to be used by an light-triggering
% trigger-box

% Close the curtains
pause(1); % Let them breathe a sec
set(annoStr.Ready,'Visible','off');

for i = 1:numTrial
    set(annoStr.plus, 'Visible', 'on');
    pause(1.0)
    set(annoStr.plus, 'Visible', 'off');
    
    set([annoStr.EEE annoStr.visTrig],'Visible','on');
    pause(1.0)
    set([annoStr.EEE annoStr.visTrig],'Visible','off');
    pause(1.0)
end
end

function endoParadigm(annoStr, numTrial)
% This sub-function mimcs the ordering of visual feedback used
% for experiments involving endoscopy recordings, and might need more
% time between trials.

FBInstr = {'Your Own Voice'; 'Loud White Noise'};

set(annoStr.curSessNote, 'Visible', 'off')
pause(1); %Let them breathe a sec
set(annoStr.Ready,'Visible','off');

for i = 1:numTrial
    curTrial = ['Trial ' num2str(i)];
    AudFB    = FBInstr(mod(i,2)+ 1); % Alternating trials    
    
    set(annoStr.trialT, 'String', curTrial);
    set(annoStr.FBCue, 'String', AudFB);
    set([annoStr.trialT annoStr.FBCue], 'Visible', 'on');
    pause() % Hit spacebar
    set([annoStr.trialT annoStr.FBCue], 'Visible', 'off');
    
    set(annoStr.plus, 'Visible', 'on');
    pause(1.0)
    set(annoStr.plus, 'Visible', 'off');
    
    set(annoStr.EEE, 'Visible','on');
    pause(1.0)
    set(annoStr.EEE, 'Visible','off');
    pause(1.0)
end
end
