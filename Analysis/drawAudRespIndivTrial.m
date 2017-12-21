function drawAudRespIndivTrial(res, plotFolder)

curSess          = res.curSess;
f0b              = round(10*res.f0b)/10;
AudFB            = res.AudFB;

time             = res.timeA;
micf0Trials      = res.audioMf0TrialPert;
heaf0Trials      = res.audioHf0TrialPert;
limits           = res.limitsA;
numTrial         = res.numPertTrialsPP;
trigs            = res.pertTrigPP;

plotpos = [10 0];
plotdim = [1200 1050];
IndivTrialAudResp = figure('Color', [1 1 1]);
set(IndivTrialAudResp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

micColor     = 'b';
headColor    = 'r';
pertBoxC     = [0.8 0.8 0.8];
fontN        = 'Arial';
legAnnoFSize = 12;
titleFSize   = 10;
axisLSize    = 10;
lineThick    = 4;

ha = tight_subplot(5, 1, [0.02 0.02],[0.05 0.12],[0.03 0.03]);

for ii = 1:numTrial      
    axes(ha(ii))
        
    pertAx  = [trigs(ii,1), trigs(ii,2)];
    pertAy  = [200 200];
    
    area(pertAx, pertAy, -200, 'FaceColor', pertBoxC, 'EdgeColor', pertBoxC)
    hold on
    mH = plot(time, micf0Trials(:, ii), micColor, 'LineWidth', lineThick);
    hold on
    hH = plot(time, heaf0Trials(:, ii), headColor, 'LineWidth', lineThick);

    if ii == 1
        xlabel('Time (s)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold')
        ylabel('f0 (cents)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold')
    end
    
    title(['Trial ' num2str(ii)], 'FontName', fontN, 'FontSize', titleFSize, 'FontWeight', 'bold')
    axis(limits); box off

    set(gca, 'FontName', fontN,...
             'FontSize', axisLSize,...
             'FontWeight','bold')
end
legend([mH hH],{'Microphone', 'Headphones'},...
            'Position', [0.8 0.30 0.1 0.1],...
            'Box', 'off',...
            'Edgecolor', [1 1 1],...
            'FontName', fontN,...
            'FontSize', legAnnoFSize,...
            'FontWeight', 'bold');
        
sup = suptitle({curSess; ['AudFB: ' AudFB]; ['f0: ' num2str(f0b) 'Hz']});
set(sup, 'FontName', fontN,...
         'FontSize', titleFSize,...
         'FontWeight','bold')
     
plots = {'AudResp_IndiTrial'};
for i = 1:length(plots)
    plTitle = [curSess '_' plots{i} '.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end