function audioDynamics_Somato = InflationResponse(secTime, secAudioMean)
% [respVar, respVarm, respVarSD, InflaStimVar] = InflationResponse(secTime, secAudio)
% Identifies the relevant pitch contour characteristics that are important
% for deciding how a participant responded to the inflation of the balloon
% during production. iR is a structure representing the result variables
% from studying the inflation response (iR). The prefix letter denotes
% whether the variable is a index (i), a time (t), or a value (v). 
%
% secTime:  Vector of time points corresponding to the sectioned data (numSamp)
% secAudio: 3D mat of sectioned audio (numSamp x numTrial x event)
%           The 1st 3D layer are Onset Sections
%           The 2nd 3D later are Offset Sections
%
% respVar: Matrix of per trial iR results (numTrial x 4)
%          respVar(:,1) = Time of the minimum f0 value in the sec
%          respVar(:,2) = Minimum f0 value in sec (stim magnitude)
%          respVar(:,3) = Value of f0 at end of sec (response magnitude)
%          respVar(:,4) = ABS percent change of stim and response
%          magnitude (response percentage)
% respVarM:    Vector of mean trial values from respVarm (1x4)
% respVarSD:   Vector of standard deviation of the trial values from respVar (1x4)
% InflaSimVar: Values of the mean time at stim magnitude and the mean stim magnitude

[numSamp, ~] = size(secAudioMean); % Size of the data we are dealing with

ir = initInflationResponseStruct(); % Initialize the structure that handles the variable calculations
ir.time     = secTime;              % Time Interval for the sectioned trials (-0.5->1.0s)
ir.onset    = secAudioMean(:, 1);   % f0 Trace sectioned around pert Onset.

ir.iAtOnset = find(ir.time == 0);
ir.tAtOnset = 0;                     % duh
ir.vAtOnset = ir.onset(ir.iAtOnset); % f0 value at t = 0

ir.iPostOnsetR = find(0 < ir.time & .20 >= ir.time); % Range of indices between t > 0ms and t =< 200ms;
[minOn, minIdx] = min(ir.onset(ir.iPostOnsetR));     % Minimum f0 val within PostOnsetR

% StimMag
ir.iAtMin  = ir.iPostOnsetR(minIdx);       % Indice of the min f0 value following trigger
ir.tAtMin  = ir.time(ir.iAtMin);           % Time at min f0 value following trigger
ir.vAtMin  = minOn;                        % Min f0 value in PostOnsetR
ir.stimMag = abs(ir.vAtMin - ir.vAtOnset); % Distance traveled from onset to min value

% RespMag
ir.tAtRespRange = [0.8 1.0];              % Time Values in the period 800ms to 1000ms after perturbation onset
ir.iAtRespRange = find(ir.time >= ir.tAtRespRange(1) & ir.time<= ir.tAtRespRange(2));                % Last index in section
ir.vAtRespRange = ir.onset(ir.iAtRespRange);    % Time Value when participant 'fully responded' (1.0s)
ir.tAtResp      = mean(ir.tAtRespRange);
ir.vAtResp      = mean(ir.vAtRespRange);
ir.respMag = ir.vAtResp - ir.vAtMin; % Distance traveled from min f0 value to response f0 value

% RespPer
ir.respPer = 100*(ir.respMag/ir.stimMag); % Percent change from stimMag to respMag 

% Add to the audioDynamics struct
respVarM = [ir.tAtMin ir.stimMag ir.respMag ir.respPer];
audioDynamics_Somato.respVarM = respVarM;
drawInflationResultMetrics(ir, 1, 0, 0); % Generates useful manuscript Fig
end

function ir = initInflationResponseStruct()

ir.time     = [];
ir.onset    = [];

ir.iAtOnset = []; % Index where t = 0
ir.tAtOnset = []; % Time at t = 0
ir.vAtOnset = []; % f0 value at t = 0

ir.iPostOnsetR = []; % Range of indices between t = 0ms and t = 200ms;
ir.iAtMin      = []; % Index at min f0 value in PostOnsetR
ir.tAtMin      = []; % Time at min f0 value in PostOnsetR
ir.vAtMin      = []; % Min f0 value in PostOnsetR
ir.stimMag     = []; % ir.vAtMin - ir.vAtOnset ..in a perfect world vAtOnset = 0

ir.iAtResp = []; % Index of f0 value when participant 'fully' responded...right now = last value in section
ir.tAtResp = []; % Time at f0 value when participant 'fully' responded
ir.vAtResp = []; % f0 value when participant 'fully' responded 
ir.respMag = []; % vAtResp - vAtMin   ...distance traveled
ir.respPer = []; % Percent change from stimMag to respMag
end