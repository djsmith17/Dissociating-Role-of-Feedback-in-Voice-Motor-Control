function drawJNDResults(JNDa, saveResultsDir)
% drawJNDResults() displays the results from a single subject's JND runs.
%
% This function calls the following helper functions
% -tight_subplot
% -export_fig

participant = JNDa.resJNDs(1).participant;
numRuns     = length(JNDa.resJNDs);

f0        = JNDa.resJNDs(1).f0;
JNDScoreM = JNDa.JNDScoreMean;

selectOpt = JNDa.resJNDs(1).selectOpt;

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
    % Set up some specifics of this recording
    resJND    = JNDa.resJNDs(ii);
    curRun    = resJND.run;
    trlsComp  = resJND.trialsCompleted;
    JNDScore  = resJND.JNDScore;
    catchAccu = resJND.LastSetAccuracy;
    
    distProg = resJND.distProgression;
    revTrl   = resJND.trialsAtReversals;
    revDst   = resJND.distAtReversals;    
    opt1CTrl = resJND.trialsAtCorrectOpt1;
    opt1CDst = resJND.distAtCorrectOpt1;
    opt1ITrl = resJND.trialsAtIncorrectOpt1;
    opt1IDst = resJND.distAtIncorrectOpt1;
    opt2CTrl = resJND.trialsAtCorrectOpt2;
    opt2CDst = resJND.distAtCorrectOpt2;
    opt2ITrl = resJND.trialsAtIncorrectOpt2;
    opt2IDst = resJND.distAtIncorrectOpt2;
    
    scoreNote = ['JND Score: ' num2str(JNDScore) ' cents'];
    accurNote = ['Last Trials Accuracy: ' num2str(catchAccu) '%'];
    
    axes(ha(ii))
    rV = plot(revTrl, revDst,'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],'MarkerSize',12);
    hold on;
    plot(distProg, 'Color', [0.7 0.7 0.7], 'LineWidth', 3);
    
    c1 = plot(opt1CTrl, opt1CDst,'o','MarkerFaceColor',tColors(1,:),'MarkerEdgeColor',tColors(1,:),'MarkerSize',10);
    i1 = plot(opt1ITrl, opt1IDst,'o','MarkerFaceColor',tColors(2,:),'MarkerEdgeColor',tColors(2,:),'MarkerSize',10);
    c2 = plot(opt2CTrl, opt2CDst,'o','MarkerFaceColor',tColors(3,:),'MarkerEdgeColor',tColors(3,:),'MarkerSize',10);
    i2 = plot(opt2ITrl, opt2IDst,'o','MarkerFaceColor',tColors(4,:),'MarkerEdgeColor',tColors(4,:),'MarkerSize',10);
       
    jS = line([0 trlsComp], [JNDScore JNDScore],'LineStyle', '-.', 'LineWidth',3,'color',[1 0 1]);
    title(curRun, 'FontSize', titleFS, 'FontName', 'Arial', 'FontWeight', 'bold')
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
    
    if ~isempty(c1)
       aH(1) = c1;
    end
    if ~isempty(i1)
       aH(2) = i1;
    end
    if ~isempty(rV)
       aH(3) = rV;
    end
    if ~isempty(c2)
       aH(4) = c2;
    end
    if ~isempty(i2)
       aH(5) = i2;
    end  
end

suptitle({'f0 Acuity JND', [participant ', f0 = ' num2str(f0) ' Hz'], ['JND Score: ' num2str(JNDScoreM) ' Cents']})

legend(aH,{['Correct ' selectOpt{1}],['Incorrect ' selectOpt{1}],'Reversals',['Correct ' selectOpt{2}],['Incorrect ' selectOpt{2}]},...
          'Orientation','Horizontal',...
          'FontSize', 12,...
          'Position', [0.51 0.455 0.02 0.05],...
          'EdgeColor', [0.5 0.5 0.5])

plTitle = [participant 'JNDStaircaseResults.jpg'];
saveFileName = fullfile(saveResultsDir, plTitle);
export_fig(saveFileName) 
end