function drawIntralTrialT(x, y, fs, span)
h = figure();

figpos = [300 200];
figdim = [1000 500];

set(h,'color',[1 1 1],'position',[figpos figdim])
t = 0:1/fs:((length(x)-1)/fs);
dottedStartx = span(1)/fs*ones(1,50);
dottedEndx   = span(2)/fs*ones(1,50);
dottedy = 5*linspace(-1,1,50);

subplot(1,2,1)
plot(t,x)
hold on
plot(dottedStartx, dottedy,'.k')
hold on
plot(dottedEndx, dottedy,'.k')
xlabel('Time (s)','FontSize', 10, 'FontWeight',  'bold')
ylabel('Amplitude (V)','FontSize', 10, 'FontWeight',  'bold')
title('Microphone', 'FontSize', 10)
axis([0 4 -0.5 0.5])
box off
set(gca, 'FontSize', 10,...
        'FontWeight', 'bold')

subplot(1,2,2)
plot(t,y)
hold on
% plot(dottedStartx, dottedy, '.k')
hold on
% plot(dottedEndx, dottedy,'.k')
xlabel('Time (s)','FontSize', 10, 'FontWeight',  'bold')
ylabel('Amplitude (V)','FontSize', 10, 'FontWeight',  'bold')
title('Headphone', 'FontSize', 10)
axis([0 3 -0.5 0.5])
box off
set(gca, 'FontSize', 10,...
        'FontWeight', 'bold')

% suptitle(['Aspiration Noise Gain of ' num2str(pert)])
end