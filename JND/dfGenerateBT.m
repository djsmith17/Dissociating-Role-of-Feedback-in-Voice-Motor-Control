function [bTRS, fs] = dfGenerateBT(dirs)
%Extract and create a baseline token for use in a JND experiment. Token is
%created from existing voice data (previous experiments)

SavFileDir    = dirs.SavFileDir;
baseTokenFile = dirs.baseTokenFile;

load(SavFileDir);
thisData  = DRF.rawData(9);      % Take the 9th trial. It will be a control trial
fsRec     = DRF.expParam.sRateAnal;
fs        = 44100;
sample    = thisData.signalIn;   % Grab the microphone channel.
  
auto        = 1;
tokenL      = .5; %stimulus duration to be played in the JND
tokenLP     = tokenL*fsRec;
riseTime    = .05;
fallTime    = .05;
riseTimeP   = riseTime*fsRec;
fallTimeP   = fallTime*fsRec;
riseQperiod = (4*riseTime)^-1;
fallQperiod = (4*fallTime)^-1;
sampLen     = length(sample);
recDuration = sampLen/fsRec;
time        = linspace(0, recDuration, sampLen);

window = ones(1, tokenLP);
window(1:riseTimeP) = sin(2*pi*riseQperiod*linspace(0, riseTime, riseTimeP)).^2;
window(tokenLP-fallTimeP + 1:tokenLP) = cos(2*pi*fallQperiod*linspace(0, fallTime, fallTimeP)).^2;

if auto == 1
    stT = 2.0;
    ix1 = fsRec*stT;
    ix2 = ix1 + tokenLP - 1;
else
    figure
    plot(time, sample, 'b'); ylim([-1 1])
    
    [x, y] = ginput(1);
    ix1 = round(x(1)*fsRec); %Choose a single point on the line with roughly .5s following it
    ix2 = ix1 + tokenLP - 1;
end

bT  = sample(ix1:ix2); %baseline Token
bTW = bT.*window';     %windowed

%resample to 44.1kHz for simplification
bTRS = resample(bTW, fs, fsRec); %resampled

%See line 71 of GAP-Pitch....
audiowrite(baseTokenFile, bTRS, fs)
end