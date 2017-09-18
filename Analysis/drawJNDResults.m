function drawJNDResults(JNDa, dirs, runs2Analyze, allRunData, allMeanJND, allCatchAcc)

saveResultsDir = dirs.SavResultsDir;

plotpos = [420 263];
plotdim = [840 525];
xyFS    = 12;
titleFS = 12;

annoPos = [.29 .77;
           .77 .77;
           .29 .31;
           .77 .31]; 

AllJND = figure('Color', [1 1 1]);
set(AllJND, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

ha = tight_subplot(2,2,[0.15 0.1],[0.1 0.1],[0.1 0.05]);

for ii = runs2Analyze
    UD = allRunData(ii);
    meanJND = allMeanJND(ii);
    catchAccu = allCatchAcc(ii);
    if isfield(UD, 'JNDDirection')
        note = UD.JNDDirection;
        anno2 = ['Catch Accuracy: ' num2str(catchAccu) '%'];
    else
        note = [];
        anno2 = ['Trials Performed: ' num2str(UD.performedTrials) '/' num2str(UD.totalTrials)];
    end
    
    axes(ha(ii))
    plot(UD.x, 'LineWidth', 3);
    hold on;
    plot(find(UD.response==1),UD.x(find(UD.response==1)),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'MarkerSize',10);
    plot(find(UD.response==0),UD.x(find(UD.response==0)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'MarkerSize',10);
    line([0 length(UD.response)], [meanJND meanJND],'LineStyle', '-.', 'LineWidth',3,'color',[1 0 1])
    title(['f0 Acuity JND ' num2str(ii) ': ' note], 'FontSize', titleFS, 'FontName', 'Arial', 'FontWeight', 'bold')
    xlabel('Trials','FontSize', xyFS,'FontName','Arial', 'FontWeight', 'bold');
    ylabel('f0 Distance (cents)','FontSize', xyFS,'FontName','Arial', 'FontWeight', 'bold');
    
    t = annotation('textbox',[annoPos(ii,1) annoPos(ii,2) 0.45 0.1],...
                   'string', {['Mean JND Score: ' num2str(meanJND) ' cents'];...
                                                  anno2},...
                    'LineStyle','none',...
                    'FontWeight','bold',...
                    'FontSize',8,...
                    'FontName','Arial');
    
    set(gca, 'LineWidth',2, 'FontSize',12,...
              'FontName','Arial',...
              'FontWeight','bold')
end

suptitle(JNDa.participant)

plTitle = [JNDa.participant 'JNDStaircaseResults.jpg'];     
saveFileName = fullfile(saveResultsDir, plTitle);
export_fig(saveFileName) 
end