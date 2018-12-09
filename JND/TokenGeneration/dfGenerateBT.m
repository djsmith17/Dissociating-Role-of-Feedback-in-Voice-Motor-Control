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

baseTokenFile = dirs.baseTokenFile;

BVDataMic = GT.BVDataMic;
BVfs      = GT.BVfs;

fs        = GT.fs;         % The sampling rate to play the tokens at.
tokenLenP = GT.tokenLenP;  % The token length in points.

% Resample the microphone recording to the sampling rate we intend to play
% it at
BVDataMicRS = resample(BVDataMic, fs, BVfs);

recLen  = length(BVDataMicRS);     % points
recDur  = recLen/fs;           % seconds
time    = linspace(0, recDur, recLen); % used for manual plotting

riseTime    = .05;               % seconds
fallTime    = .05;               % seconds
riseTimeP   = riseTime*fs;    % points
fallTimeP   = fallTime*fs;    % points
riseQperiod = (4*riseTime)^-1;
fallQperiod = (4*fallTime)^-1;

window = ones(1, tokenLenP);
window(1:riseTimeP) = sin(2*pi*riseQperiod*linspace(0, riseTime, riseTimeP)).^2;
window(tokenLenP-fallTimeP + 1:tokenLenP) = cos(2*pi*fallQperiod*linspace(0, fallTime, fallTimeP)).^2;

if auto == 1 % Automatic sectioning
    stT = 2.0;
    ix1 = fs*stT;
    ix2 = ix1 + tokenLenP - 1;
else         % Manual sectioning
    figure
    plot(time, BVDataMicRS, 'b'); ylim([-1 1])
    
    [x, ~] = ginput(1);
    ix1 = round(x(1)*fs); % Choose a single point with roughly .5s following it
    ix2 = ix1 + tokenLenP - 1;
    close all
end

BTraw = BVDataMicRS(ix1:ix2);    % raw baseline token (BT)
BT    = BTraw.*window';          % windowed 

audiowrite(baseTokenFile, BT, fs) % saved to tokens folder 
end