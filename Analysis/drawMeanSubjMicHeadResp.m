function drawMeanSubjMicHeadResp(poolRes, targPixDim, plotFolder, fLabel)

curSess          = poolRes.curSess;

numPerturb       = poolRes.numPertTrialsFin;

time             = poolRes.secTime;
micf0            = poolRes.audioMf0MeanPert{1};
headf0           = poolRes.audioMf0MeanCont;

limits  = poolRes.limitsAmean;
pltName = poolRes.pltName;

legLines = [];
legNames = {};

% Plotting Variables
plotpos        = [10 100];
plotdim        = targPixDim;
MeanSubjf0Resp = figure('Color', [1 1 1]);
set(MeanSubjf0Resp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

dottedStartx = [0 0];
dottedy      = [-500 500];
fontN        = 'Arial';
legAnnoFSize = 25;
titleFSize   = 35;
axisLSize    = 30;
lineThick    = 4;

ha = tight_subplot(1,2,[0.1 0.03],[0.12 0.15],[0.05 0.05]);

%Onset of Perturbation
axes(ha(1))

plot(dottedStartx, dottedy,'color',[0.3 0.3 0.3],'LineWidth',lineThick)
hold on

nH = shadedErrorBar(time, headf0(:,1), headf0(:,2), 'lineprops', 'r', 'transparent', 1);
set(nH.mainLine, 'LineWidth', lineThick)
hold on

nM = shadedErrorBar(time, micf0(:,1), micf0(:,2), 'lineprops', 'b', 'transparent', 1);
set(nM.mainLine, 'LineWidth', lineThick)
hold on

xlabel('Time (s)',   'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold'); 
ylabel('f0 (cents)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold')
title('Onset of Perturbation', 'FontName', fontN, 'FontSize', titleFSize, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold')

%Offset of Perturbation
axes(ha(2))

plot(dottedStartx, dottedy,'color',[0.3 0.3 0.3],'LineWidth',lineThick)
hold on

fH = shadedErrorBar(time, headf0(:,3), headf0(:,4), 'lineprops', 'r', 'transparent', 1);
set(fH.mainLine, 'LineWidth', lineThick)
hold on

fM = shadedErrorBar(time, micf0(:,3), micf0(:,4), 'lineprops', 'b', 'transparent', 1);
set(fM.mainLine, 'LineWidth', lineThick)
hold on

xlabel('Time (s)',   'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold'); 
ylabel('f0 (cents)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold')
title('Offset of Perturbation', 'FontName', fontN, 'FontSize', titleFSize, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold',...
        'YAxisLocation', 'right');

sup = suptitle(curSess);
set(sup, 'FontName', fontN,...
         'FontSize', titleFSize,...
         'FontWeight','bold')

if fLabel == 1
    figureL = pltName(end);
    figureMark = annotation('textbox', [0.01 0.88 0.05 0.1],...
                            'string', figureL,...
                            'LineStyle', 'none',...
                            'FontName', fontN,...
                            'FontSize', titleFSize,...
                            'FontWeight','bold');
end

legend(legLines, legNames,...
       'Position', [0.83 0.75 0.1 0.1],...
       'Box', 'off',...
       'Edgecolor', [1 1 1],...
       'FontName', fontN,...
       'FontSize', legAnnoFSize,...
       'FontWeight', 'bold');

plots = {'Figure'};
for i = 1:length(plots)
    plTitle = [pltName 'MicvHead.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end

function anno = checkSig(stat, thresh, anno)

if stat < thresh
    anno = [anno '*'];
end
end