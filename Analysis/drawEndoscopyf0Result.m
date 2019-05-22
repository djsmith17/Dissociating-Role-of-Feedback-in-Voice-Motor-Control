function drawEndoscopyf0Result()

close all
dirs = dfDirs('Dissociating-Role-of-Feedback-in-Voice-Motor-Control');

allParti = {'DRF5', 'DRF9', 'DRF12', 'DRF14', 'DRF19'};
eachTrial = [3 10 8 5 8];

for ii = 1:5
    participant = allParti{ii};
    run         = 'SFL1';
    curTrial    = eachTrial(ii);
    analyzeAndDrawResult(dirs, participant, run, curTrial)
end
end

function analyzeAndDrawResult(dirs, participant, run, curTrial)
        
resultParticipantDir = fullfile(dirs.Results, participant, run);
behavResultsFile     = fullfile(resultParticipantDir, [participant run 'ResultsDRF.mat']);
videoCodResultsFile  = fullfile(resultParticipantDir, [participant run 'EndoFrameMeasuresDJS.mat']);

load(behavResultsFile) % returns res
load(videoCodResultsFile) % returns CodedEndoFrameDataSet

time             = res.timef0;
sigs             = res.audioMf0TrialPert;
trialNums        = res.allIdxFin(res.pertIdxFin);
pertTrig         = res.pertTrigsFin;
limits           = res.limitsA;

timePres         = res.timeS;
sensorP          = res.sensorPsv;
pressureLim      = res.limitsP;

pertColor = [0.8 0.8 0.8];

ii = find(trialNums == curTrial);

plotpos = [10 0];
plotdim = [1000 500];
InterTrialf0 = figure('Color', [1 1 1]);
set(InterTrialf0, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

pertAx  = [pertTrig(ii,1), pertTrig(ii,2)];
pertAy  = [600 600];

% f0 Results
subplot(2,1,1)
yyaxis right
plot(timePres, sensorP(:,ii), '--k', 'LineWidth', 1.5)
ylabel('Pressure (psi)')
axis(pressureLim);
set(gca,'FontSize', 14,...
        'FontWeight','bold')
yyaxis left

pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'EdgeColor', pertColor);
hold on    

plot(time, sigs(:,ii), 'b', 'LineWidth', 2)
ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')
title({[participant ' ' run],[ ' Trial ' num2str(trialNums(ii))]}, 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off
lgd = legend(pA, 'Perturbation Period');
lgd.EdgeColor = 'none';

set(gca,'FontSize', 14,...
        'FontWeight','bold')
  
% Video Results:
curTable = CodedEndoFrameDataSet{ii};

timeFrames = linspace(0,4,121);
pt1X = curTable.FidPt1X;
pt2X = curTable.FidPt2X;
pt1Y = curTable.FidPt1Y;
pt2Y = curTable.FidPt2Y;
distRaw = curTable.Dist;
[dist, distLim] = adjustDistMeasure(timeFrames, distRaw);

XPtsDiff = pt2X - pt1X;
YPtsDiff = pt2Y - pt1Y;

% %X Component
% subplot(2,1,2)
% pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'EdgeColor', pertColor);
% hold on
% pt1XLn = plot(timeFrames, pt1X, 'm');
% hold on
% pt2XLn = plot(timeFrames, pt2X, 'g');
% axis([0 4 180 480]); box off
% ylabel('X Position (Pixel)')
% title('X Component over time')
% legend([pt1XLn, pt2XLn], 'Point1', 'Point2', 'box', 'off', 'Orientation', 'horizontal')
% 
% set(gca,'FontSize', 14,...
%         'FontWeight','bold')
% 
% subplot(2,1,3)
% pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'EdgeColor', pertColor);
% hold on
% pt1YLn = plot(timeFrames, pt1Y, 'm');
% hold on
% pt2YLn = plot(timeFrames, pt2Y, 'g');
% axis([0 4 150 350]); box off
% set(gca,'Ydir','reverse')
% ylabel('Y Position (Pixel)')
% title('Y Component over time')
% legend([pt1YLn, pt2YLn], 'Point1', 'Point2', 'box', 'off', 'Orientation', 'horizontal')
% 
% set(gca,'FontSize', 14,...
%         'FontWeight','bold')

subplot(2,1,2)
pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'EdgeColor', pertColor);
hold on
plot([-1 5], [0 0], '--k')
ptDistLn = plot(timeFrames, dist, 'b');
axis([0 4 distLim]); box off
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold')
ylabel('Distance (Pixels)')
title('Euclidian Distance Between Points')

set(gca,'FontSize', 14,...
        'FontWeight','bold')

fileName = 'EndoDistanceMeasure.png';
savedFigFile = fullfile(resultParticipantDir, [participant run fileName]);
export_fig(savedFigFile)
end

function [dist, distLim] = adjustDistMeasure(time, distRaw)

indPrePert = time < 1.0;

distSmooth = smooth(distRaw, 3);
distBase   = mean(distSmooth(indPrePert));

dist = distSmooth - distBase;
distMax = max(dist) + 10;
distMin = min(dist) - 10;
distLim = [distMin distMax];
end