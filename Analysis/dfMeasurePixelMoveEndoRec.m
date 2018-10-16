function dfMeasurePixelMoveEndoRec()
matlab.video.read.UseHardwareAcceleration('off')

close all
enA.project     = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
enA.participant = 'DRF_ENP3';    % List of multiple participants.
enA.run         = 'SFL1';
enA.ext         = 'All';
enA.trial       = 'Trial1';

dirs           = dfDirs(enA.project);

dirs.parsedVideoDir  = fullfile(dirs.SavData, enA.participant, 'parsedVideo');
dirs.parsedVideoFile = fullfile(dirs.parsedVideoDir, [enA.participant 'parsedVideo' enA.trial '.avi']);
dirs.expResultsFile = fullfile(dirs.Results, enA.participant, enA.run, [enA.participant enA.run enA.ext 'ResultsDRF.mat']);

if isfile(dirs.parsedVideoFile)
    rawVObj = VideoReader(dirs.parsedVideoFile);
else
    error('Could not find Parsed Video')
end

if isfile(dirs.expResultsFile)
    load(dirs.expResultsFile) % Returns res
else
    error('Could not find Results')
end

enA.sRateExp    = res.sRate;
enA.allIdxFin   = res.allIdxFin;
enA.numTrials   = length(enA.allIdxFin);
enA.expAudioM   = res.audioMinc;
enA.expAudioH   = res.audioHinc;
enA.trialLenP   = length(enA.expAudioM);
enA.setA        = (1:enA.trialLenP) - 1;
enA.timef0      = res.timef0;
enA.expAudioMf0 = res.audioMf0TrialPert;

% Setup Video Frame Structure
[rawVStr, enA]            = setupVideoFrames(rawVObj, enA);
[enA.vidFig, enA.vidAxes] = setViewWindowParams(enA);

frameXs = [];
frameYs = [];
for ii = 1:2
    fprintf('Frame %f\n', ii)
    image(enA.vidAxes, rawVStr(ii).cdata)
    hold(enA.vidAxes, 'on')

    [x_fid1, y_fid1] = ginput(1);
    plot(enA.vidAxes, x_fid1, y_fid1, 'b*')
    [x_fid2, y_fid2] = ginput(1);
    plot(enA.vidAxes, x_fid2, y_fid2, 'b*')
    [x_fid3, y_fid3] = ginput(1);
    plot(enA.vidAxes, x_fid3, y_fid3, 'b*')
    [x_targ, y_targ] = ginput(1);
    plot(enA.vidAxes, x_targ, y_targ, 'r*')
    hold(enA.vidAxes, 'off')

    frameXs = cat(1, frameXs, [x_fid1 x_fid2 x_fid3 x_targ]);
    frameYs = cat(1, frameYs, [y_fid1 y_fid2 y_fid3 y_targ]);
    
    if ii > 1
        diffX = diff(frameXs(ii-1:ii, :));
        diffY = diff(frameYs(ii-1:ii, :));
        
        xE1 = diffX(2) - diffX(1);
        xE2 = diffX(3) - diffX(1);
        xE3 = diffX(2) - diffX(3);
        
        yE1 = diffY(2) - diffY(1);
        yE2 = diffY(1) - diffY(3);
        yE3 = diffY(2) - diffY(3);
    end
end

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

function [vidFig, vidAxes] = setViewWindowParams(enA)

vidFig  = figure();
vidAxes = axes();

set(vidFig,'position',[150 150 enA.vidWidth enA.vidHeight]);
set(vidAxes,'units','pixels');
set(vidAxes,'position',[0 0 enA.vidWidth enA.vidHeight]);
end