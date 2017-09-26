function drawJNDResults(JNDa, dirs, numRuns, allRunData, allMeanJND, allCatchAcc)

saveResultsDir = dirs.SavResultsDir;

plotpos = [50 150];
plotdim = [1600 750];
xyFS    = 12;
titleFS = 12;

tColors = [[0 0 1]; %Correct Different
           [1 0 0]; %Incorrect Different
           [0.5 0.5 1]; %Correct Same
           [1 0.5 0.5]];%Incorrect Same

annoPos = [.38 .79;
           .87 .79;
           .38 .32;
           .87 .32]; 
       
aH = [0 0 0 0 0];

AllJND = figure('Color', [1 1 1]);
set(AllJND, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

ha = tight_subplot(2,2,[0.15 0.06],[0.1 0.1],[0.05 0.03]);

for ii = 1:numRuns
    UD = allRunData(ii);
    meanJND = allMeanJND(ii);
    catchAccu = allCatchAcc(ii);
    revNote   = [num2str(UD.reversals) ' Reversals, '];
    triNote   = [num2str(UD.performedTrials) '/' num2str(UD.totalTrials) ' Trials, '];
    timNote   = [num2str(round(10*UD.elapsedTime)/10) ' min'];
    
    anno2 = ['Catch Accuracy: ' num2str(catchAccu) '%'];
    
    axes(ha(ii))
    rV = plot(find(UD.reversal~=0), UD.x(find(UD.reversal~=0)),'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],'MarkerSize',12);
    hold on;
    plot(UD.x, 'Color', [0.7 0.7 0.7], 'LineWidth', 3);   
    
    cD = plot(find(UD.allTrialTypes==1), UD.x(find(UD.allTrialTypes==1)),'o','MarkerFaceColor',tColors(1,:),'MarkerEdgeColor',tColors(1,:),'MarkerSize',10);
    iD = plot(find(UD.allTrialTypes==2), UD.x(find(UD.allTrialTypes==2)),'o','MarkerFaceColor',tColors(2,:),'MarkerEdgeColor',tColors(2,:),'MarkerSize',10);
    cS = plot(find(UD.allTrialTypes==3), UD.x(find(UD.allTrialTypes==3)),'o','MarkerFaceColor',tColors(3,:),'MarkerEdgeColor',tColors(3,:),'MarkerSize',10);
    iS = plot(find(UD.allTrialTypes==4), UD.x(find(UD.allTrialTypes==4)),'o','MarkerFaceColor',tColors(4,:),'MarkerEdgeColor',tColors(4,:),'MarkerSize',10);
       
    aJ = line([0 length(UD.response)], [meanJND meanJND],'LineStyle', '-.', 'LineWidth',3,'color',[1 0 1]);
    title(['f0 Acuity JND ' num2str(ii) ': ' revNote triNote timNote], 'FontSize', titleFS, 'FontName', 'Arial', 'FontWeight', 'bold')
    xlabel('Trials','FontSize', xyFS,'FontName','Arial', 'FontWeight', 'bold');
    ylabel('f0 Distance (cents)','FontSize', xyFS,'FontName','Arial', 'FontWeight', 'bold');
    axis([0 80 0 70])
    
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
    
    if ~isempty(cD)
       aH(1) = cD;
    end
    if ~isempty(iD)
       aH(2) = iD;
    end
    if ~isempty(rV)
       aH(3) = rV;
    end
    if ~isempty(cS)
       aH(4) = cS;
    end
    if ~isempty(iS)
       aH(5) = iS;
    end  
end

suptitle({JNDa.participant; ['f0: ' num2str(UD.subjf0) ' Hz']})

legend(aH,{['Correct ' UD.tN{1}],['Incorrect ' UD.tN{1}],'Reversals',['Correct ' UD.tN{2}],['Incorrect ' UD.tN{2}]},...
       'Orientation','Horizontal',...
       'FontSize', 8,...
       'Position', [0.51 0.48 0.01 0.01],...
       'EdgeColor', [0 0 0])


plTitle = [JNDa.participant 'JNDStaircaseResults.jpg'];     
saveFileName = fullfile(saveResultsDir, plTitle);
export_fig(saveFileName) 
end