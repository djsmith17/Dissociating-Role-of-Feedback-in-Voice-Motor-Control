function dfMeasurePixelMoveEndoRec()
matlab.video.read.UseHardwareAcceleration('off')

close all
enA.project     = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
enA.participant = 'DRF_ENP3';    % List of multiple participants.
enA.run         = 'SFL1';
enA.ext         = 'All';
enA.trial       = 1;

dirs           = dfDirs(enA.project);

dirs.parsedVideoDir  = fullfile(dirs.SavData, enA.participant, 'parsedVideo');
dirs.parsedVideoFile = fullfile(dirs.parsedVideoDir, [enA.participant 'parsedVideoTrial' num2str(enA.trial) '.avi']);
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
enA.trialLenT   = enA.trialLenP/enA.sRateExp;
enA.setA        = (1:enA.trialLenP) - 1;
enA.timef0      = res.timef0;
enA.expAudioMf0 = res.audioMf0TrialPert;
enA.pertTrigs   = res.pertTrigsFin(enA.trial, :);

videoWidth  = 720;
videoHeight = 480;
frameRate   = 30;

videoFileReader = vision.VideoFileReader(dirs.parsedVideoFile);
videoPlayer     = vision.VideoPlayer('Position',[200,500,videoWidth + 30, videoHeight+ 20]);
objectFrame     = videoFileReader();

% figure; imshow(objectFrame);
% objectRegion = round(getPosition(imrect)); close;
objectRegion = [126, 66, 457, 336];

iniPoints    = detectMinEigenFeatures(rgb2gray(objectFrame), 'ROI', objectRegion);

pointImage = insertShape(objectFrame, 'Rectangle', objectRegion, 'Color','red');
pointImage = insertMarker(pointImage, iniPoints.Location,'+','Color','white');
figure; imshow(pointImage);
specROI = round(getPosition(imrect));
specInd = (iniPoints.Location(:,1) > specROI(1) & iniPoints.Location(:,1) < specROI(1) + specROI(3)) & (iniPoints.Location(:,2) > specROI(2) & iniPoints.Location(:,2) < specROI(2) + specROI(4));

% multiStabInd = [];
% for ss = 1:3
%     stabROI = round(getPosition(imrect)); 
%     stabInd = find((iniPoints.Location(:,1) > stabROI(1) & iniPoints.Location(:,1) < stabROI(1) + stabROI(3)) & (iniPoints.Location(:,2) > stabROI(2) & iniPoints.Location(:,2) < stabROI(2) + stabROI(4)));
%     multiStabInd = cat(1, multiStabInd, stabInd);
% end
stabInd = ~specInd;
close;

tracker = vision.PointTracker('MaxBidirectionalError',1);
initialize(tracker,iniPoints.Location,objectFrame);

quivFig = figure();
quivAx  = axes(quivFig);
set(quivAx, 'XLim', [0 720], 'YLim', [0 480])
set(quivAx, 'YDir', 'reverse')

frame1        = videoFileReader();
movMean       = frame1;
correctedMean = frame1;
curFrame      = frame1;
[F1Points, ~] = tracker(curFrame);

curFrameCorrected   = frame1;
lastPointsCorrected = F1Points;

allFrames = struct;
allFrames(1).frame = frame1;
allPoints = [];
allPoints = cat(1, allPoints, F1Points);

ii = 2;
Hcumulative = eye(3);
while ~isDone(videoFileReader)
    % Read in new frame
    frame = videoFileReader();

    lastFrame = curFrame; % z^-1
    lastFrameCorrected = curFrameCorrected;
    curFrame  = frame;
    movMean   = movMean + curFrame;
    
    [curPoints, validity] = tracker(curFrame);
    [lastPoints, ~]       = tracker(lastFrame);
    
%     curPointsStab  = curPoints(stabInd, :);
%     lastPointsStab = lastPoints(stabInd, :);
    
    % Estimate transform from frame A to frame B, and fit as an s-R-t
%     H = cvexEstStabilizationTform(lastFrame, curFrame, lastPointsStab, curPointsStab);
%     HsRt = cvexTformToSRT(H);
%     Hcumulative = HsRt * Hcumulative;
%     curFrameCorrected = imwarp(curFrame,affine2d(Hcumulative),'OutputView',imref2d(size(curFrame)));
% 
%     % Detect Velocity of Points between this frame and last frame
%     [curPointsCorrected, ~]  = tracker(curFrameCorrected);
%     [lastPointsCorrected, ~] = tracker(lastFrameCorrected);
    
    svPoints  = curPoints;
    svLPoints = lastPoints;
    
    pointsShift = svPoints - svLPoints;   
    
    % Generate Quiver plot of change in points
    quiver(quivAx, curPoints(:,1), curPoints(:,2), pointsShift(:,1), pointsShift(:,2), 0.5);
    set(quivAx, 'XLim', [0 720], 'YLim', [0 480])
    set(quivAx, 'YDir', 'reverse')
    title(['Frame ' num2str(ii)])
    
    allPoints = cat(3, allPoints, svPoints);

    % Display as color composite with last corrected frame
%     videoPlayer(imfuse(lastFrameCorrected, curFrameCorrected,'ColorChannels','red-cyan'));
%     correctedMean = correctedMean + curFrameCorrected;
    
    out = insertMarker(frame, curPoints(validity, :), '+');
    videoPlayer(out);
    
    allFrames(ii).frame = frame;
    ii = ii+1;
end
nFrames = ii - 1;
% correctedMean = correctedMean/(ii-2);
% movMean = movMean/(ii-2);

release(videoPlayer);
release(videoFileReader);

% figure; imshowpair(movMean, correctedMean, 'montage');
% title(['Raw input mean', repmat(' ',[1 50]), 'Corrected sequence mean']);

videoTime = linspace(0, enA.trialLenT, nFrames);

velAllX = diff(squeeze(allPoints(:, 1, :))');
velAllY = diff(squeeze(allPoints(:, 2, :))');
velSpecX = diff(squeeze(allPoints(specInd, 1, :))');
velSpecY = diff(squeeze(allPoints(specInd, 2, :))');
velStabX = diff(squeeze(allPoints(stabInd, 1, :))');
velStabY = diff(squeeze(allPoints(stabInd, 2, :))');

figure
subplot(1,2,1)
plot(videoTime(2:end), velStabX, 'b')
hold on
plot(videoTime(2:end), velSpecX, 'r')
plot([enA.pertTrigs(1) enA.pertTrigs(1)], [-50 50], 'k--')
plot([enA.pertTrigs(2) enA.pertTrigs(2)], [-50 50], 'k--')
axis([0 4 -40 40])

subplot(1,2,2)
plot(videoTime(2:end), velStabY, 'b')
hold on
plot(videoTime(2:end), velSpecY, 'r')
plot([enA.pertTrigs(1) enA.pertTrigs(1)], [-50 50], 'k--')
plot([enA.pertTrigs(2) enA.pertTrigs(2)], [-50 50], 'k--')
axis([0 4 -40 40])
 
% figure
% subplot(1,2,1)
% plot(videoTime(2:end), velAllX)
% hold on
% plot([enA.pertTrigs(1) enA.pertTrigs(1)], [-50 50], 'k--')
% plot([enA.pertTrigs(2) enA.pertTrigs(2)], [-50 50], 'k--')
% axis([0 4 -40 40])
% 
% subplot(1,2,2)
% plot(videoTime(2:end), velAllY)
% hold on
% plot([enA.pertTrigs(1) enA.pertTrigs(1)], [-50 50], 'k--')
% plot([enA.pertTrigs(2) enA.pertTrigs(2)], [-50 50], 'k--')
% axis([0 4 -40 40])
end