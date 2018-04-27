function drawMeanSubjf0Resp(poolRes, targPixDim, plotFolder, fLabel)

curSess          = poolRes.curSess;
numControl       = poolRes.numContTrialsFin;
numPerturb       = poolRes.numPertTrialsFin;

time             = poolRes.secTime;
contf0           = poolRes.audioMf0MeanCont;  
pertF0           = poolRes.audioMf0MeanPert;

meanf0ContOnset  = contf0(:,1);
CIf0ContOnset    = contf0(:,2);
meanf0ContOffset = contf0(:,3);
CIf0ContOffset   = contf0(:,4);

meanf0PertOnsetM  = pertF0{1}(:,1);
CIf0PertOnsetM    = pertF0{1}(:,2);
meanf0PertOffsetM = pertF0{1}(:,3);
CIf0PertOffsetM   = pertF0{1}(:,4);

meanf0PertOnsetV  = pertF0{2}(:,1);
CIf0PertOnsetV    = pertF0{2}(:,2);
meanf0PertOffsetV = pertF0{2}(:,3);
CIf0PertOffsetV   = pertF0{2}(:,4);

statLib = poolRes.statLib;
statSMM = round(statLib(1), 1);
statSMV = round(statLib(2), 1);
statRMM = round(statLib(3), 1);
statRMV = round(statLib(4), 1);
statRPM = round(statLib(5));
statRPV = round(statLib(6));
statSP  = statLib(7);
statRP  = statLib(8);
statPP  = statLib(9);

limits = poolRes.limitsAmean;
pltName = poolRes.pltName;

pValueThresh = 0.05;

% Plotting Variables
plotpos        = [10 100];
plotdim        = targPixDim;
MeanSubjf0Resp = figure('Color', [1 1 1]);
set(MeanSubjf0Resp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

dottedStartx = [0 0];
dottedy      = [-500 500];
maskColor    = 'b';
voicColor    = 'r';
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
nC = shadedErrorBar(time, meanf0ContOnset, CIf0ContOnset, 'lineprops', 'k', 'transparent', 1); %Voice
hold on
nM = shadedErrorBar(time, meanf0PertOnsetM, CIf0PertOnsetM, 'lineprops', maskColor, 'transparent', 1); %Masked
hold on
nV = shadedErrorBar(time, meanf0PertOnsetV, CIf0PertOnsetV, 'lineprops', voicColor, 'transparent', 1); %Voice

set(nM.mainLine, 'LineWidth', lineThick)
set(nV.mainLine, 'LineWidth', lineThick)
set(nC.mainLine, 'LineWidth', lineThick)
xlabel('Time (s)',   'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold'); 
ylabel('f0 (cents)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold')
title('Onset of Perturbation', 'FontName', fontN, 'FontSize', titleFSize, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'XTickLabel',{'-0.5' '0' '0.5' '1.0'},...
        'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold')

%Offset of Perturbation
axes(ha(2))
plot(dottedStartx, dottedy,'color',[0.3 0.3 0.3],'LineWidth',lineThick)
hold on
fC = shadedErrorBar(time, meanf0ContOffset, CIf0ContOffset, 'lineprops', 'k', 'transparent', 1); %Voice
hold on
fM = shadedErrorBar(time, meanf0PertOffsetM, CIf0PertOffsetM, 'lineprops', maskColor, 'transparent', 1); %Masked
hold on
fV = shadedErrorBar(time, meanf0PertOffsetV, CIf0PertOffsetV, 'lineprops', voicColor, 'transparent', 1); %Voice

set(fM.mainLine, 'LineWidth', lineThick)
set(fV.mainLine, 'LineWidth', lineThick)
set(fC.mainLine, 'LineWidth', lineThick)
xlabel('Time (s)',   'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold'); 
ylabel('f0 (cents)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold')
title('Offset of Perturbation', 'FontName', fontN, 'FontSize', titleFSize, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'XTickLabel', {'-0.5' '0' '0.5' '1.0'},...
        'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold',...
        'YAxisLocation', 'right');

sup = suptitle(curSess);
set(sup, 'FontName', fontN,...
         'FontSize', titleFSize,...
         'FontWeight','bold')

% Done plotting, now to add some annotations
annoStim = ['SM (M/NM): ' num2str(statSMM) ' cents / ' num2str(statSMV) ' cents'];
annoResp = ['RM (M/NM): ' num2str(statRMM) ' cents / ' num2str(statRMV) ' cents'];
annoPerc = ['RP (M/NM): ' num2str(statRPM) '% / ' num2str(statRPV) '%'];

annoStim = checkSig(statSP, pValueThresh, annoStim);
annoResp = checkSig(statRP, pValueThresh, annoResp);
annoPerc = checkSig(statPP, pValueThresh, annoPerc);
 
statBox = annotation('textbox',[.30 .75 0.45 0.1],...
                     'string', {annoStim;
                                annoResp
                                annoPerc},...
                      'LineStyle','none',...
                      'FontName', fontN,...
                      'FontSize', legAnnoFSize,...
                      'FontWeight','bold');

if fLabel == 1
    figureL = pltName(end);
    figureMark = annotation('textbox', [0.01 0.88 0.05 0.1],...
                            'string', figureL,...
                            'LineStyle', 'none',...
                            'FontName', fontN,...
                            'FontSize', titleFSize,...
                            'FontWeight','bold');
end

legend([fC.mainLine fM.mainLine fV.mainLine],{[num2str(numControl) ' Control Trials'], [num2str(numPerturb(1)) ' Masked Trials'], [num2str(numPerturb(2)) ' Not Masked Trials']},...
            'Position', [0.83 0.75 0.1 0.1],...
            'Box', 'off',...
            'Edgecolor', [1 1 1],...
            'FontName', fontN,...
            'FontSize', legAnnoFSize,...
            'FontWeight', 'bold');

plots = {'Figure'};
for i = 1:length(plots)
    plTitle = [pltName '.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end

function anno = checkSig(stat, thresh, anno)

if stat < thresh
    anno = [anno '*'];
end
end