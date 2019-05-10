
time             = res.timef0;
sigs             = res.audioMf0TrialPert;
trialNums        = res.allIdxFin(res.pertIdxFin);
pertTrig         = res.pertTrigsFin;
limits           = res.limitsA;

timePres         = res.timeS;
sensorP          = res.sensorPsv;
pressureLim      = res.limitsP;

pertColor = [0.8 0.8 0.8];

ii = 10;

plotpos = [10 100];
plotdim = [1000 250];
InterTrialf0 = figure('Color', [1 1 1]);
set(InterTrialf0, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

pertAx  = [pertTrig(ii,1), pertTrig(ii,2)];
pertAy  = [600 600];

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
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold')
ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')
title(['Trial ' num2str(trialNums(ii))], 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'FontSize', 14,...
        'FontWeight','bold')
    

timeFrames = linspace(0,4,121);
pt1X = curTable.FidPt1X;
pt2X = curTable.FidPt2X;
pt1Y = curTable.FidPt1Y;
pt2Y = curTable.FidPt2Y;
XPtsDiff = pt2X - pt1X;
YPtsDiff = pt2Y - pt1Y;

ChangeX = figure('Color', [1 1 1]);
set(ChangeX, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'EdgeColor', pertColor);
hold on
plot(timeFrames, XPtsDiff, 'k');
axis([0 4 0 100]); box off
xlabel('Time (s)')
ylabel('Diff X (Pt2X - Pt1X)')
title('Difference in X Components')

set(gca,'FontSize', 14,...
        'FontWeight','bold')
    
ChangeY = figure('Color', [1 1 1]);
set(ChangeY, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'EdgeColor', pertColor);
hold on
plot(timeFrames, YPtsDiff, 'k');
axis([0 4 -250 -160]); box off
xlabel('Time (s)')
ylabel('Diff Y (Pt2Y - Pt1Y)')
title('Difference in Y Components')

set(gca,'FontSize', 14,...
        'FontWeight','bold')