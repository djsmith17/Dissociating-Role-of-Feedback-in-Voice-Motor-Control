function expParam = dfInitExpParam()

% Top level experiment
expParam.project      = [];
expParam.expType      = [];
expParam.curSess      = [];
expParam.rmsB         = [];

% The Subject
expParam.subject      = [];
expParam.gender       = []; % Has default
expParam.age          = []; % Has default
expParam.f0b          = []; % Has default
expParam.targRMS      = []; % Has default

% The Run
expParam.run          = [];
expParam.balloon      = [];
expParam.tightness    = [];
expParam.InflaVarNm   = [];
expParam.elapsedTime  = [];
expParam.loudResults  = []; 

% NIDAQ Info
expParam.niDev        = []; % NIDAQ Device Name. For more information, see dfInitNIDAQ
expParam.sRateQ       = []; % NIDAQ sampling rate
expParam.niCh         = []; % Structure of Channel Names

% Audapter Info
expParam.sRate              = [];
expParam.frameLen           = [];
expParam.downFact           = [];
expParam.sRateAnal          = [];
expParam.frameLenDown       = [];
expParam.audioInterfaceName = [];
expParam.ostFN              = [];
expParam.pcfFN              = [];

% Experimental Specifics
expParam.trialLen     = []; % Seconds
expParam.numTrial     = []; % Number of trials
expParam.curTrial     = [];
expParam.perCatch     = []; % Percentage of catch or 'perturbed' trials
expParam.numMaskRep   = [];

% Auditory Feedback Settings
expParam.headGain     = []; % Output gain above the input
expParam.AudFB        = [];
expParam.AudFBSw      = [];
expParam.AudPert      = [];
expParam.AudPertSw    = [];

% Experimental Timing Settings
expParam.rdyPause  = [];
expParam.cuePause  = [];
expParam.buffPause = [];
expParam.endPause  = [];
expParam.resPause  = [];
expParam.boundsRMS = [];

% Set some default parameters
expParam = setDefaultExpParam(expParam);
end

function expParam = setDefaultExpParam(expParam)

expParam.project = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType = 'Somatosensory Perturbation_Perceptual';

expParam.gender  = 'female';
expParam.age     = 20;
expParam.f0b     = 200;
expParam.targRMS = 70;

expParam.balloon   = 'N/A';
expParam.tightness = 'N/A';

expParam.niDev   = 'Dev2';

expParam.trialLen = 4; %seconds

expParam.headGain  = 5; %dB SPL
expParam.AudFB     = 'Voice Shifted';
expParam.AudFBSw   = 1;
expParam.AudPert   = '-100 cents ramped';
expParam.AudPertSw = 1;

expParam.InflaVarNm = 'N/A';

expParam.rdyPause  = 5.0; % How long to show them 'Ready'
expParam.cuePause  = 1.0; % How long the cue period lasts
expParam.buffPause = 0.8; % Give them a moment to start speaking
expParam.endPause  = 0.5;
expParam.resPause  = 2.0; % How long the rest/VisFB lasts
expParam.boundsRMS = 3;   % +/- dB
end