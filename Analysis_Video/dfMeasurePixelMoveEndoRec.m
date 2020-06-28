function dfMeasurePixelMoveEndoRec()
matlab.video.read.UseHardwareAcceleration('off')

close all
enA.project     = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
enA.participant = 'DRF12';    % List of multiple participants.
enA.run         = 'SFL1';
enA.ext         = 'All';
enA.trial       = 6;
enA.curSess     = [enA.participant 'parsedVideoTrial' num2str(enA.trial)];
monitorSize     = get(0, 'Monitor');
enA.monSize     = monitorSize(1, 3:4);

dirs           = dfDirs(enA.project);

dirs.parsedVideoDir  = fullfile(dirs.SavDataEndo, enA.participant, 'parsedVideo');
dirs.parsedVideoFile = fullfile(dirs.parsedVideoDir, [enA.curSess '.avi']);
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

enA.videoWidth  = 720;
enA.videoHeight = 480;
enA.frameRate   = 30;

videoPlayerPos = [(enA.monSize(1)/2 - enA.videoWidth - 50), 400, (enA.videoWidth + 30), (enA.videoHeight+ 20)];
quiverPos      = [(enA.monSize(1)/2 + 50), 400, (enA.videoWidth + 30), (enA.videoHeight+ 20)];
videoFileReader = vision.VideoFileReader(dirs.parsedVideoFile);
videoPlayer     = vision.VideoPlayer('Position', videoPlayerPos);
objectFrame     = videoFileReader();

% figure; imshow(objectFrame);
% objectRegion = round(getPosition(imrect)); close;
objectRegion = [126, 66, 457, 336];

[compassLines, textPos] = setupVideoFrameAnno();

objectFrame = insertShape(objectFrame, 'Line', compassLines, 'Color', 'white');
objectFrame = insertText(objectFrame, textPos, {'Left', 'Anterior'}, 'TextColor', 'white', 'BoxColor', 'black', 'AnchorPoint', 'CenterBottom');

iniPoints    = detectMinEigenFeatures(rgb2gray(objectFrame), 'ROI', objectRegion);

pointImage = insertShape(objectFrame, 'Rectangle', objectRegion, 'Color','red');
pointImage = insertMarker(pointImage, iniPoints.Location,'+','Color','white');
figure('Position', videoPlayerPos); imshow(pointImage);
specROI = round(getPosition(imrect));
specInd = (iniPoints.Location(:,1) > specROI(1) & iniPoints.Location(:,1) < specROI(1) + specROI(3)) & (iniPoints.Location(:,2) > specROI(2) & iniPoints.Location(:,2) < specROI(2) + specROI(4));

% multiStabInd = [];
% for ss = 1:3
%     stabROI = round(getPosition(imrect)); 
%     stabInd = find((iniPoints.Location(:,1) > stabROI(1) & iniPoints.Location(:,1) < stabROI(1) + stabROI(3)) & (iniPoints.Location(:,2) > stabROI(2) & iniPoints.Location(:,2) < stabROI(2) + stabROI(4)));
%     multiStabInd = cat(1, multiStabInd, stabInd);
% end
% stabInd = multiStabInd;
stabInd = ~specInd;
close;

tracker = vision.PointTracker('MaxBidirectionalError',1);
initialize(tracker,iniPoints.Location,objectFrame);

quivFig = figure('Position', quiverPos);
quivAx  = axes(quivFig);
set(quivAx, 'XLim', [0 720], 'YLim', [0 480])
set(quivAx, 'YDir', 'reverse')

frame1        = videoFileReader();
movMean       = frame1;
correctedMean = frame1;
curFrame          = frame1;
curFrameCorrected = frame1;
[F1Points, ~] = tracker(curFrame);

allFrames = struct;
allFrames(1).frame = frame1;
allPoints = [];
allPoints = cat(1, allPoints, F1Points);

stable = 0;
% videoPlayer(frame1)
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
    
    if stable == 1
        curPointsStab  = curPoints(stabInd, :);
        lastPointsStab = lastPoints(stabInd, :);

        % Estimate transform from frame A to frame B, and fit as an s-R-t
        H = cvexEstStabilizationTform(lastFrame, curFrame, lastPointsStab, curPointsStab);
        HsRt = cvexTformToSRT(H);
        Hcumulative = HsRt * Hcumulative;
        curFrameCorrected = imwarp(curFrame,affine2d(Hcumulative),'OutputView',imref2d(size(curFrame)));

        % Detect Velocity of Points between this frame and last frame
        [curPointsCorrected, ~]  = tracker(curFrameCorrected);
        [lastPointsCorrected, ~] = tracker(lastFrameCorrected);
    end
    
    svPoints  = curPoints;
    svLPoints = lastPoints;
    
    pointsShift = svPoints - svLPoints;   
    
    % Generate Quiver plot of change in points
    quiver(quivAx, curPoints(:,1), curPoints(:,2), pointsShift(:,1), pointsShift(:,2), 0.5);
    set(quivAx, 'XLim', [0 720], 'YLim', [0 480])
    set(quivAx, 'YDir', 'reverse')
    title(['Frame ' num2str(ii)])
    
    allPoints = cat(3, allPoints, svPoints);
    
    if stable == 1
        % Display as color composite with last corrected frame
        out = imfuse(lastFrameCorrected, curFrameCorrected,'ColorChannels','red-cyan');
        correctedMean = correctedMean + curFrameCorrected;
    else
        out = frame;
    end
    
    out = insertMarker(out, curPoints(validity, :), '+');
    out = insertShape(out, 'Line', compassLines, 'Color', 'white');
    out = insertText(out, textPos, {'Left', 'Anterior'}, 'TextColor', 'white', 'BoxColor', 'black', 'AnchorPoint', 'CenterBottom');     
    
    videoPlayer(out);
    
    allFrames(ii).frame = out;
    ii = ii+1;
end
close(quivFig)
enA.nFrames = ii - 1;
correctedMean = correctedMean/(ii-2);
movMean = movMean/(ii-2);

release(videoPlayer);
release(videoFileReader);

alteredVideoFile = fullfile(dirs.Results, enA.participant, enA.run, [enA.curSess 'PixelTrackVideo.avi']);

vidWriter = VideoWriter(alteredVideoFile, 'Uncompressed AVI');
vidWriter.FrameRate = 30;
open(vidWriter);

for nn = 1:enA.nFrames
    writeVideo(vidWriter, allFrames(nn).frame)
end
close(vidWriter)

if stable == 1
    figure; imshowpair(movMean, correctedMean, 'montage');
    title(['Raw input mean', repmat(' ',[1 50]), 'Corrected sequence mean']);
end
plotPointVelocities(dirs, enA, allPoints, specInd, stabInd)
end

function [compassLines, textPos] = setupVideoFrameAnno()

compassSt = [640 380];
vertArrow = [compassSt(1) compassSt(2) compassSt(1) compassSt(2)+50];
horzArrow = [compassSt(1) compassSt(2) compassSt(1)+50 compassSt(2)];
vertArrowUp = [vertArrow(3) vertArrow(4) vertArrow(3)-5 vertArrow(4)-5];
vertArrowDn = [vertArrow(3) vertArrow(4) vertArrow(3)+5 vertArrow(4)-5];
horzArrowUp = [horzArrow(3) horzArrow(4) horzArrow(3)-5 horzArrow(4)+5];
horzArrowDn = [horzArrow(3) horzArrow(4) horzArrow(3)-5 horzArrow(4)-5];

compassLines = [vertArrow; horzArrow; vertArrowUp; vertArrowDn; horzArrowUp; horzArrowDn];

leftTextPos = [mean([horzArrow(1) horzArrow(3)]), (horzArrow(2) -2)];
AntTextPos  = [vertArrow(1), (vertArrow(4) + 25)];
textPos     = [leftTextPos; AntTextPos];
end

function plotPointVelocities(dirs, enA, allPoints, specInd, stabInd)
videoTime = linspace(0, enA.trialLenT, enA.nFrames);

velAllX = diff(squeeze(allPoints(:, 1, :))')./enA.frameRate;
velAllY = diff(squeeze(allPoints(:, 2, :))')./enA.frameRate;
velSpecX = velAllX(:, specInd);
velSpecY = velAllY(:, specInd);
velStabX = velAllX(:, stabInd);
velStabY = velAllY(:, stabInd);

velXMin = min(min(velAllX)); velYMin = min(min(velAllY));
velXMax = max(max(velAllX)); velYMax = max(max(velAllY));

limitsX = [0 4 velXMin-0.1 velXMax + 0.1];
limitsY = [0 4 velYMin-0.1 velYMax + 0.1];

figSize = [1400 500];
figPos  = [(enA.monSize(1)/2-figSize(1)/2) 90];

sepGroupsFig = figure('Color', [1 1 1]);
set(sepGroupsFig, 'Position', [figPos figSize])

ha = tight_subplot(1,2,[0.08 0.08],[0.12 0.15],[0.05 0.01]);

axes(ha(1))
plot(videoTime(2:end), velStabX, 'b')
hold on
plot(videoTime(2:end), velSpecX, 'r')
plot([enA.pertTrigs(1) enA.pertTrigs(1)], [-50 50], 'k--')
plot([enA.pertTrigs(2) enA.pertTrigs(2)], [-50 50], 'k--')

x = [0.1 0.1]; y = [0.7 0.85];
annotation('textarrow',x,y,'String','Left')
axis(limitsX); box off
xlabel('Time (s)')
ylabel('Velocity (pixels/s)')
title('Horizontal Component of Pixel Velocity')

axes(ha(2))
plot(videoTime(2:end), velStabY, 'b')
hold on
plot(videoTime(2:end), velSpecY, 'r')
plot([enA.pertTrigs(1) enA.pertTrigs(1)], [-50 50], 'k--')
plot([enA.pertTrigs(2) enA.pertTrigs(2)], [-50 50], 'k--')

x = [0.62 0.62]; y = [0.7 0.85];
annotation('textarrow',x,y,'String','Anterior')
axis(limitsY); box off
xlabel('Time (s)')
ylabel('Velocity (pixels/s)')
title('Vertical Component of Pixel Velocity')

enA.curSess(strfind(enA.curSess, '_')) = '';
suptitle(enA.curSess)

saveFileName = fullfile(dirs.Results, enA.participant, enA.run, [enA.curSess '_PixelVelocity.jpg']);
export_fig(saveFileName)
end