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

%Paradigm Configurations
expParam.project    = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType    = 'Somatosensory Perturbation_Perceptual';
expParam.gender     = 'male';
expParam.trialLen   = 4;  % Seconds
expParam.numTrial   = 1;
expParam.AudFBSw    = 2;    % Masking Noise

% Set our dirs based on the project
dirs = dfDirs(expParam.project);

% Paradigm Configurations
expParam.sRate              = 48000;  % Hardware sampling rate (before downsampling)
expParam.frameLen           = 96;     % Number of samples in processing frame (before downsampling)
expParam.downFact           = 3;
expParam.sRateAnal          = expParam.sRate/expParam.downFact;
expParam.frameLenDown       = expParam.frameLen/expParam.downFact;
expParam.audioInterfaceName = 'MOTU MicroBook'; %'ASIO4ALL' 'Komplete'

%Set up Audapter
Audapter('deviceName', expParam.audioInterfaceName);
Audapter('setParam', 'downFact', expParam.downFact, 0);
Audapter('setParam', 'sRate', expParam.sRateAnal, 0);
Audapter('setParam', 'frameLen', expParam.frameLenDown, 0);
p = getAudapterDefaultParams(expParam.gender);

%Set up OST and PCF Files.
expParam.ostFN = fullfile(dirs.Prelim, 'SFPerturbOST.ost'); check_file(expParam.ostFN);
expParam.pcfFN = fullfile(dirs.Prelim, 'SFPerturbPCF.pcf'); check_file(expParam.pcfFN);

% Set up Auditory Feedback (Masking Noise)
[expParam, p]      = dfSetAudFB(expParam, dirs, p);

% Pause Periods
expParam.readyPause = 5.0;
expParam.cuePause  = 1.0; % How long the cue period lasts
expParam.buffPause = 0.8; %Give them a moment to start speaking
expParam.endPause  = 0.5;
expParam.resPause  = 2.0; % How long the rest/VisFB lasts

noiseLen = calcMaskLen(expParam);

ET = tic;
%Open the curtains
pause(expParam.readyPause);                % Let them breathe a sec

for ii = 1:expParam.numTrial
    
    %Set the OST and PCF functions
    Audapter('ost', expParam.ostFN, 0);
    Audapter('pcf', expParam.pcfFN, 0);
    
    pause(expParam.cuePause)
    
    fprintf('Trial %d\n', ii)
    AudapterIO('init', p);
    Audapter('reset');
    Audapter('start');
    pause(expParam.buffPause)

    pause(expParam.trialLen)
    pause(expParam.endPause)
    
    Audapter('stop');
    pause(expParam.resPause)
end
elapsed_time = toc(ET);   % Elapsed Time of the experiment
fprintf('\nElapsed Time: %f (s)\n', elapsed_time)
end

function noiseLen = calcMaskLen(expParam)

numTrial = expParam.numTrial;

rdyTime = expParam.readyPause;
cueTime = expParam.cuePause;
buffTime = expParam.buffPause;
trlTime  = expParam.trialLen;
endTime  = expParam.endPause;
resTime  = expParam.resPause;

noiseLen = rdyTime + (cueTime + buffTime + trlTime + endTime + resTime)*numTrial;
end