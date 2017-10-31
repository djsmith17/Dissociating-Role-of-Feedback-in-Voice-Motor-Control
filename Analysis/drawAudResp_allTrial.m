function drawAudResp_AllTrial(res, curSess, curRec, plotFolder)

time  = res.time;
trigs = res.allTrialTrigs;
trialf0_mic  = res.allTrialf0(:,1,:);
trialf0_head = res.allTrialf0(:,2,:);
baselinef0 = res.meanTrialf0b;
numTrial   = res.trialCount;
limits     = res.f0Limits;

numRow = numTrial/5;

plotpos = [100 100];
plotdim = [1800 800];
InterTrialAudResp = figure('Color', [1 1 1]);
set(InterTrialAudResp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

curSess(strfind(curSess, '_')) = ' ';
curRec(strfind(curRec, '_')) = ' ';

pertColor = [0.8 0.8 0.8];

ha = tight_subplot(numRow, 5, [0.1 0.02],[0.05 0.12],[0.05 0.03]);

for ii = 1:numTrial      
    axes(ha(ii))
        
    pertAx  = [trigs(ii,1), trigs(ii,2)];
    pertAy  = [200 200];
    
    area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor)
    hold on
    tM = plot(time, trialf0_mic(:, ii), 'b', 'LineWidth', 1.5);
    hold on
    tH = plot(time, trialf0_head(:, ii), 'r', 'LineWidth', 1.5);

    if ii == 1
        xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')
    end
    
    title(['Trial ' num2str(ii)], 'FontSize', 10, 'FontWeight', 'bold')
    axis(limits); box off

    set(gca,'FontSize', 10,...
            'FontWeight','bold')
end
l0 = legend([tM tH], 'Microphone', 'Headphones'); 
set(l0,'box', 'off','FontSize', 14, 'FontWeight', 'bold');
        

suptitle({curSess; [curRec '   f0: ' num2str(baselinef0) 'Hz']})

plots = {'AllTrialf0_AudResp'};
for i = 1:length(plots)
    plTitle = [curSess '_' plots{i} '_' curRec '.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end