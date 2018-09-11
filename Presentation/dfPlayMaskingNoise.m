function dfPlayMaskingNoise()
%Test the levels of masking noise based on how the paradigm is set up for
%voice perturbations
%
% This script calls the following 8 functions:
% dfDirs.m
% dfSetAudFB.m
%
% This uses the toolbox from MATLAB-Toolboxes
% speechres

recordedRMS = (5+rand(1991,1))/1000;
rmsMean     = mean(recordedRMS);

expParam = dfInitExpParam();

expParam.curSess = 'Masking Noise Test';

expParam.numTrial   = 1;
expParam.numMaskRep = expParam.numTrial;
expParam.AudFB      = 'AC Masking Noise';
expParam.AudFBSw    = 2; % Masking Noise

% Set our dirs based on the project
dirs = dfDirs(expParam.project);

% Set up Audapter Parameters and function paths
[expParam, p] = dfInitAudapater(dirs, expParam);

% Set up Auditory Feedback (Masking Noise)
[p, SSNw, SSNfs] = dfSetAudFB(expParam, dirs, p);

% This is where the fun begins
fprintf('\nStarting Trials\n\n')

% Dim the lights (Set the visual Feedback)
[msrStr, annoStr] = dfSetVisFB(1, expParam.curSess, expParam.targRMS, expParam.boundsRMS);

% Only play masking noise for this condition
if expParam.AudFBSw == 2
    sound(SSNw, SSNfs)
end

ET = tic;
% Open the curtains
pause(expParam.rdyPause); % Let them breathe a sec
set(annoStr.Ready, 'Visible','off'); % Turn off 'Ready?'
for ii = 1:expParam.numTrial
    
    %Set the OST and PCF functions
    Audapter('ost', expParam.ostFN, 0);
    Audapter('pcf', expParam.pcfFN, 0);
    
    %Cue to begin trial
    set(annoStr.plus, 'Visible','on');
    pause(expParam.cuePause)
    set(annoStr.plus, 'Visible','off');
    
    %Phonation Start
    set([annoStr.EEE annoStr.visTrig],'Visible','on');
    
    fprintf('Trial %d\n', ii)
    AudapterIO('init', p);
    Audapter('reset');
    Audapter('start');
    pause(expParam.buffPause)

    pause(expParam.trialLen)
    
    %Phonation End
    set([annoStr.EEE annoStr.visTrig],'Visible','off');
    pause(expParam.endPause)
    Audapter('stop');
    
    % Compare against baseline and updated Visual Feedback
    [color, newPos, loudResult] = dfUpdateVisFB(msrStr, rmsMean);
     
    % Provide Loudness Feedback
    set(annoStr.LoudRec, 'position', newPos);
    set(annoStr.LoudRec, 'Color', color, 'FaceColor', color);
    set([annoStr.LoudRec annoStr.fbLines], 'Visible', 'on');    
    pause(expParam.resPause)
    set([annoStr.LoudRec annoStr.fbLines], 'Visible', 'off');
end
close all

elapsed_time = toc(ET);   % Elapsed Time of the experiment
fprintf('\nElapsed Time: %f (s)\n', elapsed_time)
end