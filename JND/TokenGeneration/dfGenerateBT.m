function BT = dfGenerateBT(dirs, GT, varargin)
% dfGenerateBT loads a baseline voice recording and creates a baseline 
% voice token to be used by a JND experiment.
%
% INPUTS:
% dirs:      Struc containing where the baseline voice is located and where
%            we should save the baseline token.
% GT:        Structure of variables related to generating speech tokens.
% varargin:  Flag to turn on/off automated baseline voice file sectioning.
%            1(default) is automated. 0 is manual.
%
% OUTPUTS:
% BT:        Baseline token signal

if isempty(varargin)
    auto = 1;
else
    auto = varargin{1};
end

baseTrial = GT.baseTrial; % Trial number from the baseline recordings to use.
tokenL    = GT.tokenLen;  % The target token length
fs        = GT.fs;        % The sampling rate we want the tokens played at.

baseTokenFile = dirs.baseTokenFile; % Where to save the .wav file

if ~exist(dirs.BaseFile, 'file')
    error('ERROR: No baseline voice file at the designated location!')
else
    fprintf('Loading Baseline recording from...\n%s\n\n', dirs.BaseFile)
    load(dirs.BaseFile) % Returns a structure called 'DRF'
end
    
BVData  = DRF.rawData(baseTrial); % Grab a specific trial
baseRec = BVData.signalIn;        % Grab the microphone channel.
fsRec   = DRF.expParam.sRateAnal; % Rate it was recorded at
recLen  = length(baseRec);        % points
recDur  = recLen/fsRec;           % seconds
time    = linspace(0, recDur, recLen); % used for manual plotting

tokenLP     = tokenL*fsRec;      % points at rate of recording, not audio play rate. 
riseTime    = .05;               % seconds
fallTime    = .05;               % seconds
riseTimeP   = riseTime*fsRec;    % points
fallTimeP   = fallTime*fsRec;    % points
riseQperiod = (4*riseTime)^-1;
fallQperiod = (4*fallTime)^-1;

window = ones(1, tokenLP);
window(1:riseTimeP) = sin(2*pi*riseQperiod*linspace(0, riseTime, riseTimeP)).^2;
window(tokenLP-fallTimeP + 1:tokenLP) = cos(2*pi*fallQperiod*linspace(0, fallTime, fallTimeP)).^2;

if auto == 1 % Automatic sectioning
    stT = 2.0;
    ix1 = fsRec*stT;
    ix2 = ix1 + tokenLP - 1;
else         % Manual sectioning
    figure
    plot(time, baseRec, 'b'); ylim([-1 1])
    
    [x, ~] = ginput(1);
    ix1 = round(x(1)*fsRec); % Choose a single point with roughly .5s following it
    ix2 = ix1 + tokenLP - 1;
    close all
end

BTraw = baseRec(ix1:ix2);         % raw baseline token (BT)
BTw   = BTraw.*window';           % windowed 
BT    = resample(BTw, fs, fsRec); % resampled at 44.1kHz

audiowrite(baseTokenFile, BT, fs) % saved to tokens folder 
end