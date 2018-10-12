function dfParseRawVideoScript()

close all
enA.project     = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
enA.participant = 'DRF_ENP3';    % List of multiple participants.
enA.run         = 'SFL1';
enA.ext         = 'All';

dirs           = dfDirs(enA.project);

dirs.rawVideoFile = fullfile(dirs.SavData, enA.participant, 'rawVideo', [enA.participant ' rawVideo.avi']);  % Where to find data
dirs.rawAudioFile = fullfile(dirs.SavData, enA.participant, 'rawVideo', [enA.participant ' rawAudio.wav']);  % Where to find data
dirs.parsedVideoDir = fullfile(dirs.SavData, enA.participant, 'parsedVideo');

if ~exist(dirs.parsedVideoDir, 'dir')
    mkdir(dirs.parsedVideoDir)
end

dirs.expResultsFile = fullfile(dirs.Results, enA.participant, enA.run, [enA.participant enA.run enA.ext 'ResultsDRF.mat']);

if isfile(dirs.rawVideoFile)
    disp('Loading...')
    rawVObj = VideoReader(dirs.rawVideoFile);
    [endoAudio, endoAfs] = audioread(dirs.rawAudioFile);
else
    error('Not a raw video file')
end

if isfile(dirs.expResultsFile)
    load(dirs.expResultsFile) % Returns res
else
    error('Could not find Results')
end

enA.sRateExp    = res.sRateExp;
enA.allIdxFin   = res.allIdxFin;
enA.numTrials   = length(enA.allIdxFin);
enA.expAudioM   = res.rawAudioM(:,enA.allIdxFin);
enA.expAudioH   = res.rawAudioH(:,enA.allIdxFin);
enA.trialLenP   = length(enA.expAudioM);
enA.setA        = (1:enA.trialLenP) - 1;

% Setup Video Frame Structure
[rawVStr, enA]            = setupVideoFrames(rawVObj, enA);
[enA.vidFig, enA.vidAxes] = setViewWindowParams(enA);

enA.timef0      = res.timef0;
enA.expAudioMf0 = res.audioMf0TrialPert;

enA.endoAudioDN = resample(endoAudio(:,1), enA.sRateExp, endoAfs);

enA.trialDelays     = zeros(enA.numTrials, 1);
enA.trialRecStTimes = zeros(enA.numTrials, 1);
enA.videoFrameSeg   = zeros(enA.numTrials, 2);
enA.endoAudioM      = zeros(enA.trialLenP, enA.numTrials);
for ii = 1:enA.numTrials
    enA.trialDelays(ii)     = xCorrTimeLag(enA.endoAudioDN, enA.expAudioM(:,ii));
    enA.trialRecStTimes(ii) = round(enA.trialDelays(ii)/enA.sRateExp, 2);
    enA.endoAudioM(:,ii)    = enA.endoAudioDN(enA.setA + enA.trialDelays(ii));
    
    videoFrameSt = floor(enA.trialRecStTimes(ii)*enA.vidFrameR);
    videoFrameSp = videoFrameSt + enA.trialLenF;
    
    enA.videoFrameSeg(ii, :) = [videoFrameSt videoFrameSp];
    
    trialVideo = rawVStr(videoFrameSt:videoFrameSp);
%     playRawVideo(enA, trialVideo)
%     figure
%     subplot(2,1,1)
%     plot(enA.expAudioM(:,ii))
%     subplot(2,1,2)
%     plot(enA.endoAudioM(:,ii))
%     title([num2str(enA.trialRecStTimes(ii)) 's'])
    curVideoFile = fullfile(dirs.parsedVideoDir, [enA.participant 'parsedVideoTrial' num2str(enA.allIdxFin(ii)) '.avi']);
    curV = VideoWriter(curVideoFile);
    open(curV)
    writeVideo(curV, trialVideo);
    close(curV)
end
parseTrialT = table(enA.trialRecStTimes, enA.videoFrameSeg(:, 1), 'VariableNames',{'Time', 'Frame'});

dirs.rawVideoParseNoteFile = fullfile(dirs.SavData, enA.participant, 'rawVideo', [enA.participant ' parseNotes.txt']);
writetable(parseTrialT, dirs.rawVideoParseNoteFile)
% fullVid = rawVStr;
% playRawVideo(enA, fullVid)
end

function enA = initEndoAnalysisStr()



end

function [rawVStr, enA] = setupVideoFrames(rawVObj, enA)

enA.vidWidth  = rawVObj.Width;
enA.vidHeight = rawVObj.Height;
enA.vidFrameR = rawVObj.FrameRate;
enA.trialLenF = round(enA.trialLenP*(enA.vidFrameR/enA.sRateExp));
enA.setF      = (1:enA.trialLenF) - 1;

rawVStr = struct('cdata', zeros(enA.vidHeight, enA.vidWidth, 3, 'uint8'), 'colormap', []);

enA.nFrames = 0;
while hasFrame(rawVObj)
    enA.nFrames = enA.nFrames+1;
    rawVStr(enA.nFrames).cdata = readFrame(rawVObj);
end
end

function pointLag = xCorrTimeLag(sig1, sig2)
% xCorrTimeLag(sig1, sig2, fs) calculates the lag between two (seemingly) 
% identical time based signals. 
%
% if timeLag is positive, then sig1 leads sig2. 
% if timeLag is negative, then sig1 lags sig2.

% Simple crosscorrelation between two signals
% Finds the largest peak of the result
[r, lags]    = xcorr(sig1, sig2);
[~, peakInd] = max(r);
pointLag     = lags(peakInd);
end

function [vidFig, vidAxes] = setViewWindowParams(enA)

vidFig  = figure();
vidAxes = axes();

set(vidFig,'position',[150 150 enA.vidWidth enA.vidHeight]);
set(vidAxes,'units','pixels');
set(vidAxes,'position',[0 0 enA.vidWidth enA.vidHeight]);
end

function playRawVideo(enA, vid)

movie(enA.vidAxes, vid, 1, enA.vidFrameR)
end

