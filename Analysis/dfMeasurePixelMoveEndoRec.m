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

if ~isfile(dirs.parsedVideoFile)
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

videoFileReader = vision.VideoFileReader(dirs.parsedVideoFile);
videoPlayer     = vision.VideoPlayer('Position',[200,500,750,500]);
objectFrame     = videoFileReader();

figure; imshow(objectFrame);
objectRegion = round(getPosition(imrect)); close;
iniPoints    = detectMinEigenFeatures(rgb2gray(objectFrame), 'ROI', objectRegion);

pointImage = insertMarker(objectFrame, iniPoints.Location,'+','Color','white');
figure; imshow(pointImage);
specROI = round(getPosition(imrect));
specInd = find((iniPoints.Location(:,1) > specROI(1) & iniPoints.Location(:,1) < specROI(1) + specROI(3)) & (iniPoints.Location(:,2) > specROI(2) & iniPoints.Location(:,2) < specROI(2) + specROI(4)));
stabROI = round(getPosition(imrect)); 
stabInd = (iniPoints.Location(:,1) > stabROI(1) & iniPoints.Location(:,1) < stabROI(1) + stabROI(3)) & (iniPoints.Location(:,2) > stabROI(2) & iniPoints.Location(:,2) < stabROI(2) + stabROI(4));
close;

tracker = vision.PointTracker('MaxBidirectionalError',1);
initialize(tracker,iniPoints.Location,objectFrame);

quivFig = figure();
quivAx  = axes(quivFig);
set(quivAx, 'XLim', [0 720], 'YLim', [0 480])
set(quivAx, 'YDir', 'reverse')

stabPointsIni = iniPoints.Location(stabInd,:);

i = 1;
allPoints = [];
while ~isDone(videoFileReader)
    frame = videoFileReader();
    [points, validity] = tracker(frame);    
    
    if i > 1
        stabPointsFrame = points(stabInd,:);
        stabDiff        = stabPointsFrame - stabPointsIni;
        meanDiff        = mean(stabDiff, 1);
        
        newPoints = points - meanDiff;
        
        allDiff  = newPoints - pointsIni;
        quiver(quivAx, newPoints(:,1), newPoints(:,2), allDiff(:,1), allDiff(:,2), 0.5);
        set(quivAx, 'XLim', [0 720], 'YLim', [0 480])
        set(quivAx, 'YDir', 'reverse')
        
        pointsIni     = newPoints;
        allPoints = cat(3, allPoints, newPoints);
    else
        pointsIni = points;
    end
    i = i + 1;
    
    
    out = insertMarker(frame, points(validity, :), '+');
    videoPlayer(out);    
end

release(videoPlayer);
release(videoFileReader);

% % Setup Video Frame Structure
% [rawVStr, enA]            = setupVideoFrames(rawVObj, enA);
% [enA.vidFig, enA.vidAxes] = setViewWindowParams(enA);
% 
% frameXs = [];
% frameYs = [];
% targXChange = zeros(60, 1);
% targYChange = zeros(60, 1);
% for ii = 1:60
%     fprintf('Frame %f\n', ii)
%     image(enA.vidAxes, rawVStr(ii).cdata)
%     hold(enA.vidAxes, 'on')
% 
%     [x_fid1, y_fid1] = ginput(1);
%     plot(enA.vidAxes, x_fid1, y_fid1, 'b*')
% %     [x_fid2, y_fid2] = ginput(1);
% %     plot(enA.vidAxes, x_fid2, y_fid2, 'b*')
% %     [x_fid3, y_fid3] = ginput(1);
% %     plot(enA.vidAxes, x_fid3, y_fid3, 'b*')
%     [x_targ, y_targ] = ginput(1);
%     plot(enA.vidAxes, x_targ, y_targ, 'r*')
%     hold(enA.vidAxes, 'off')
% 
%     frameXs = cat(1, frameXs, [x_fid1 x_targ]);
%     frameYs = cat(1, frameYs, [y_fid1 y_targ]);
%     
%     if ii > 1
%         diffX = diff(frameXs(ii-1:ii, :));
%         diffY = diff(frameYs(ii-1:ii, :));
%         
%         targXChange(ii) = diffX(2) + diffX(1);
%         targYChange(ii) = diffY(2) + diffX(1);
%     end
% end
% figure
% subplot(2,1,1)
% plot(targXChange)
% subplot(2,1,2)
% plot(targYChange)
% 
% figure
% plot(frameXs(:,1), frameYs(:,1), 'b*')
% hold on
% plot(frameXs(:,2), frameYs(:,2), 'r*')
% axis([0 rawVObj.Width 0 rawVObj.Height])
% set(gca, 'Ydir','reverse')
% 
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