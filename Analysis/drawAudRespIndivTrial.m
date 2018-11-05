function drawAudRespIndivTrial(res, plotFolder)

curSess          = res.curSess;
f0b              = round(res.f0b, 1); % Baseline f0 rounded to 0.1 Hz
AudFB            = res.AudFB;
numPT            = res.numPertTrialsFin;
pertTrig         = res.pertTrigsFin;
MHDelays         = res.allAuMHDelays;
AuNiDelays       = round(res.allAuNiDelays, 3);

time             = res.timef0;
micf0Trials      = res.audioMf0TrialPert;
heaf0Trials      = res.audioHf0TrialPert;
limits           = res.limitsAudRes;

plotpos = [10 10];
plotdim = [1600 800];
IndivTrialAudResp = figure('Color', [1 1 1]);
set(IndivTrialAudResp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

micColor     = 'b';
headColor    = 'r';
pertBoxC     = [0.8 0.8 0.8];
fontN        = 'Arial';
legAnnoFSize = 12;
titleFSize   = 14;
axisLSize    = 14;
lineThick    = 4;

numColums = 5;
numrows   = ceil(numPT/numColums);

ha = tight_subplot(numrows, numColums, [0.15 0.05],[0.12 0.15],[0.05 0.03]);

for ii = 1:numPT      
    axes(ha(ii))
        
    pertAx  = [pertTrig(ii,1), pertTrig(ii,2)];
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
    
    title({['Trial ' num2str(ii)], [num2str(pertTrig(ii,1)) '  ' num2str(pertTrig(ii,2))]}, 'FontName', fontN, 'FontSize', titleFSize, 'FontWeight', 'bold')
    axis(limits); box off

    set(gca, 'FontName', fontN,...
             'FontSize', axisLSize,...
             'FontWeight','bold')
end
legend([mH hH], {'Microphone', 'Headphones'},...
                 'Position', [0.8 0.93 0.05 0.05],...
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