function drawAllTrialf0(time, allTrialf0, runTrialOrder, trigs, limits, meanTrialf0b, curExp, curRecording, plotFolder)
plotpos = [50 100];
plotdim = [1500 800];
AllTrialf0 = figure('Color', [1 1 1]);
set(AllTrialf0, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

curExp(strfind(curExp, '_')) = ' ';
curRecording(strfind(curRecording, '_')) = ' ';

pertRuns = find(runTrialOrder == 1);
numPertRuns = length(pertRuns);
pertColor = [0.8 0.8 0.8];

ha = tight_subplot(5,2,[0.10 0.07],[0.10 0.10],[0.05 0.01]);

for ii = 1:numPertRuns    
    axes(ha(ii))
    
    pertAx  = [trigs(ii,1), trigs(ii,2)];
    pertAy  = [200 200];
    
    pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on
    
    fT = plot(time,allTrialf0(:,1,pertRuns(ii)), 'LineWidth',2);
    xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 12, 'FontWeight', 'bold')
    title(['Trial ' num2str(ii)], 'FontSize', 14, 'FontWeight', 'bold')
    axis([0 4 -100 100]); box off
    
    set(gca,'FontSize', 14,...
            'FontWeight','bold')
% l0 = legend([uH.mainLine pH.mainLine],[num2str(counts(1)) ' Control Trials'], [num2str(counts(2)) ' Perturb Trials']); 
% set(l0,'box', 'off','FontSize', 14, 'FontWeight', 'bold');

end

suptitle({[curExp ': Mic Recording, Perturbed Trials']; [curRecording '   f0: ' num2str(meanTrialf0b) 'Hz']})
legend([pA fT], {'Perturbation Period', 'f0 Trace'},...
                'Orientation','Horizontal',...
                'FontSize', 8,...
                'Position', [0.77 0.95 0.01 0.01],...
                'EdgeColor', [0 0 0])

plots = {'AllTrialf0'};
for i = 1:length(plots)
    plTitle = [curExp '_' plots{i} '.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end