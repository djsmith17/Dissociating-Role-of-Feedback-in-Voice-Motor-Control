function drawJNDResults(JNDa, dirs, allRunData, allMeanJND)

saveResultsDir = dirs.SavResultsDir;

plotpos = [420 263];
plotdim = [840 525];
xyFS    = 15;

AllJND = figure('Color', [1 1 1]);
set(AllJND, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

ha = tight_subplot(2,2,[0.1 0.1],[0.10 0.10],[0.1 0.1]);

for ii = 1:4
    UD = allRunData(ii);
    meanJND = allMeanJND(ii);
    
    axes(ha(ii))
    plot(UD.x, 'LineWidth', 3);
    hold on;
    plot(find(UD.response==1),UD.x(find(UD.response==1)),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'MarkerSize',10);
    plot(find(UD.response==0),UD.x(find(UD.response==0)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'MarkerSize',10);
    line([0 length(UD.response)], [meanJND meanJND],'LineStyle', '-.', 'LineWidth',3,'color',[1 0 1])
    title(['JND Acuity Run ' num2str(ii)])
    xlabel('Trials','FontSize', xyFS,'FontName','Arial');
    ylabel('Perturbations','FontSize', xyFS,'FontName','Arial');
    box off
end

suptitle(JNDa.participant)

plTitle = [JNDa.participant 'JNDStaircaseResults.jpg'];     
saveFileName = fullfile(saveResultsDir, plTitle);
export_fig(saveFileName) 
end