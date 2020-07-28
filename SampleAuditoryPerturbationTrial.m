plotPos = [-1700 500];
plotDim = [900 350];

MHFig = figure('Color', [1 1 1]);
set(MHFig, 'Position', [plotPos plotDim],'PaperPositionMode','auto')

VOTime     = 0.1;
pPertTime  = 1.0;
pertLength = 1.5;

OnsetTime = VOTime + pPertTime;
OffsetTime = OnsetTime + pertLength;

fs = 44100;
t = 0:1/fs:4;
OnsetPoint = round(OnsetTime*fs);
OffsetPoint = round(OffsetTime*fs);

signal = zeros(size(t));

rampOnPoints = round(0.11*fs);
rampOfPoints = round(0.15*fs);

rampOnPert = linspace(0, -100, rampOnPoints);
rampOfPert = linspace(-100, 0, rampOfPoints);

rampOnEnd = OnsetPoint + rampOnPoints -1;
rampOfEnd = OffsetPoint + rampOfPoints -1;

signal(OnsetPoint:rampOnEnd) = rampOnPert;
signal(OffsetPoint:rampOfEnd) = rampOfPert;
signal(rampOnEnd:OffsetPoint) = -100;

%Analysis Windows
OnsetWinAx = [-0.5 1.0]+OnsetTime;
OffsetWinAx = [-0.5 1.0]+OffsetTime;
WinAy = [600 600];

area(OnsetWinAx, WinAy, -600, 'FaceColor', 'g', 'EdgeColor', 'g', 'FaceAlpha', 0.25, 'EdgeAlpha', 0.25)
hold on
area(OffsetWinAx, WinAy, -600, 'FaceColor', 'r', 'EdgeColor', 'r', 'FaceAlpha', 0.25, 'EdgeAlpha', 0.25)

plot(t, signal, 'LineWidth', 2)
hold on
xlabel('Time (s)')
ylabel('f0 Shift (cents)')


plot([VOTime VOTime], [-600 600], '--m', 'LineWidth', 2)

plot([OnsetTime OnsetTime], [-600 600], '--g', 'LineWidth', 2)

plot([OffsetTime OffsetTime], [-600 600], '--r', 'LineWidth', 2)

box off
axis([0 4 -110 10])
set(gca, 'FontName', 'Arial',...
         'FontSize', 14,...
         'FontWeight','bold')