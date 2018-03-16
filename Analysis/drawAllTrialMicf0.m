function drawAllTrialMicf0(res, plotFolder)

curSess          = res.curSess;
AudFB            = res.AudFB;
f0b              = round(10*res.f0b)/10;
f0Type           = res.f0Type;
etMH             = res.etMH;
numPT            = res.numPertTrialsPP;
pertTrig         = res.pertTrigPP;

time             = res.timef0;
sigs             = res.audioMf0TrialPert;
limits           = res.limitsA;

curSess(strfind(curSess, '_')) = ' ';

plotpos = [10 100];
plotdim = [1600 800];
InterTrialf0 = figure('Color', [1 1 1]);
set(InterTrialf0, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

pertColor = [0.8 0.8 0.8];

ha = tight_subplot(2,5,[0.15 0.05],[0.12 0.15],[0.08 0.08]);

for ii = 1:numPT
    axes(ha(ii))
    
    pertAx  = [pertTrig(ii,1), pertTrig(ii,2)];
    pertAy  = [600 600];
    
    pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on    
    plot(time, sigs(:,ii), 'b', 'LineWidth', 2)
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold')
    ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')
    title(['Perturbed Trial ' num2str(ii)], 'FontSize', 18, 'FontWeight', 'bold')
    axis(limits); box off
    
    set(gca,'FontSize', 14,...
            'FontWeight','bold')
    
end

sup = suptitle({curSess; [ 'AudFB: ' AudFB]; ['f0: ' num2str(f0b) ' Hz']});
set(sup, 'FontSize', 18, 'FontWeight', 'bold')

timeBox = annotation('textbox',[.80 .88 0.45 0.1],...
                     'string', {f0Type;
                            ['Analysis Time: ' num2str(etMH) ' min']},...
                        'LineStyle','none',...
                        'FontWeight','bold',...
                        'FontSize',8,...
                        'FontName','Arial');

plots = {'IntraTrialf0'};
for i = 1:length(plots)
    plTitle = [curSess '_' plots{i} '.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end