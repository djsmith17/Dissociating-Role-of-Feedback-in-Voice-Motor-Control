function drawJNDResults(JNDa, saveResultsDir, allRunData)
% drawJNDResults() displays the results from a single subject's JND runs.
%
% This function calls the following helper functions
% -tight_subplot
% -export_fig

runs    = JNDa.runs;
numRuns = length(runs);

JNDScores  = JNDa.JNDScores;
lastSetAcc = JNDa.lastSetAccuracy;

plotpos = [50 150];
plotdim = [1600 750];
xyFS    = 12;
titleFS = 12;

tColors = [[0 0 1];      % Correct Different
           [1 0 0];      % Incorrect Different
           [0.5 0.5 1];  % Correct Same
           [1 0.5 0.5]]; % Incorrect Same

annoPos = [.32 .785;
           .81 .785;
           .32 .315;
           .81 .315]; 
       
aH = [0 0 0 0 0];

AllJND = figure('Color', [1 1 1]);
set(AllJND, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

ha = tight_subplot(2, 2, [0.15 0.06], [0.1 0.1], [0.05 0.03]);

for ii = 1:numRuns
    UD        = allRunData(ii);
    JNDScore  = JNDScores(ii);
    catchAccu = lastSetAcc(ii);
%     revNote   = [num2str(UD.reversals) ' Reversals, '];
%     triNote   = [num2str(UD.performedTrials) '/' num2str(UD.totalTrials) ' Trials, '];
%     timNote   = [num2str(round(10*UD.elapsedTime)/10) ' min, '];
%     
%     if isfield(UD, 'inst')
%        timNote   = [timNote 'Instru: ' UD.inst];
%     end
    
    scoreNote = ['JND Score: ' num2str(JNDScore) ' cents'];
    accurNote = ['Last Trials Accuracy: ' num2str(catchAccu) '%'];
    
    axes(ha(ii))
    rV = plot(find(UD.reversal~=0), UD.x(find(UD.reversal~=0)),'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],'MarkerSize',12);
    hold on;
    plot(UD.x, 'Color', [0.7 0.7 0.7], 'LineWidth', 3);   
    
    cD = plot(find(UD.allTrialTypes==1), UD.x(find(UD.allTrialTypes==1)),'o','MarkerFaceColor',tColors(1,:),'MarkerEdgeColor',tColors(1,:),'MarkerSize',10);
    iD = plot(find(UD.allTrialTypes==2), UD.x(find(UD.allTrialTypes==2)),'o','MarkerFaceColor',tColors(2,:),'MarkerEdgeColor',tColors(2,:),'MarkerSize',10);
    cS = plot(find(UD.allTrialTypes==3), UD.x(find(UD.allTrialTypes==3)),'o','MarkerFaceColor',tColors(3,:),'MarkerEdgeColor',tColors(3,:),'MarkerSize',10);
    iS = plot(find(UD.allTrialTypes==4), UD.x(find(UD.allTrialTypes==4)),'o','MarkerFaceColor',tColors(4,:),'MarkerEdgeColor',tColors(4,:),'MarkerSize',10);
       
    aJ = line([0 length(UD.response)], [JNDScore JNDScore],'LineStyle', '-.', 'LineWidth',3,'color',[1 0 1]);
%     title(['f0 Acuity JND ' num2str(ii) ': ' revNote triNote timNote], 'FontSize', titleFS, 'FontName', 'Arial', 'FontWeight', 'bold')
    title(runs{ii}, 'FontSize', titleFS, 'FontName', 'Arial', 'FontWeight', 'bold')
    xlabel('Trials','FontSize', xyFS,'FontName','Arial', 'FontWeight', 'bold');
    ylabel('f0 Distance (cents)','FontSize', xyFS,'FontName','Arial', 'FontWeight', 'bold');
    axis([0 80 0 80])
    
    t = annotation('textbox',[annoPos(ii,1) annoPos(ii,2) 0.45 0.1],...
                   'string', {scoreNote; accurNote},...
                    'LineStyle','none',...
                    'FontWeight','bold',...
                    'FontSize',14,...
                    'FontName','Arial');
    
    set(gca, 'LineWidth', 2,...
             'FontSize', 12,...
             'FontName', 'Arial',...
             'FontWeight', 'bold')
    
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

suptitle({'f0 Acuity JND', [JNDa.participant ', f0 = ' num2str(JNDa.f0) ' Hz'], ['JND Score: ' num2str(JNDa.JNDScoreMean) ' Cents']})

legend(aH,{['Correct ' JNDa.tN{1}],['Incorrect ' JNDa.tN{1}],'Reversals',['Correct ' JNDa.tN{2}],['Incorrect ' JNDa.tN{2}]},...
          'Orientation','Horizontal',...
          'FontSize', 12,...
          'Position', [0.51 0.455 0.02 0.05],...
          'EdgeColor', [0.5 0.5 0.5])


plTitle = [JNDa.participant 'JNDStaircaseResults' JNDa.runType '.jpg'];     
saveFileName = fullfile(saveResultsDir, plTitle);
export_fig(saveFileName) 
end