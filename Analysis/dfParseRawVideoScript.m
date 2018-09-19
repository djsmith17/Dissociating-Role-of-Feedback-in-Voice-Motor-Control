function dfParseRawVideoScript()

close all
enA.project     = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
enA.participant = 'DRF_ENP1';    % List of multiple participants.

dirs           = dfDirs(enA.project);

dirs.rawVideoFile = fullfile(dirs.SavData, enA.participant, 'rawVideo', [enA.participant ' rawVideo.avi']);  % Where to find data

if isfile(dirs.rawVideoFile)
    disp('Loading...')
    rawVObj = VideoReader(dirs.rawVideoFile);
else
    error('Not a raw video file')
end

enA.vidWidth  = rawVObj.Width;
enA.vidHeight = rawVObj.Height;
enA.vidFrameR = rawVObj.FrameRate;

rawVStr = struct('cdata', zeros(enA.vidHeight, enA.vidWidth, 3, 'uint8'), 'colormap', []);

enA.nFrames = 0;
while hasFrame(rawVObj)
    enA.nFrames = enA.nFrames+1;
    rawVStr(enA.nFrames).cdata = readFrame(rawVObj);
end

playRawVideo(rawVObj, rawVStr)
end

function playRawVideo(rawVObj, rawVStr)

vidWidth  = rawVObj.Width;
vidHeight = rawVObj.Height;
vidFrameR = rawVObj.FrameRate;

vidFig  = figure();
vidAxes = axes();

set(vidFig,'position',[150 150 vidWidth vidHeight]);
set(vidAxes,'units','pixels');
set(vidAxes,'position',[0 0 vidWidth vidHeight]);

movie(rawVStr, 1, vidFrameR)
end