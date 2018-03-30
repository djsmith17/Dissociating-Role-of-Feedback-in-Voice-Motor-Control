function drawMeanSubjf0Resp_CollarPosDiff(CR, targPixDim, pltName, plotFolder)

curSess       = CR.curSess;
numCont       = CR.numContTrials;
numPertLP     = CR.numPertTrialsLP;
numPertuLP    = CR.numPertTrialsuLP;
numPertCC     = CR.numPertTrialsCC;

time             = CR.secTime;
meanf0ContOnset  = CR.audioMf0MeanCont(:,1);
CIf0ContOnset    = CR.audioMf0MeanCont(:,2);
meanf0ContOffset = CR.audioMf0MeanCont(:,3);
CIf0ContOffset   = CR.audioMf0MeanCont(:,4);

meanf0PertOnsetLP  = CR.audioMf0MeanPertLP(:,1);
CIf0PertOnsetLP    = CR.audioMf0MeanPertLP(:,2);
meanf0PertOffsetLP = CR.audioMf0MeanPertLP(:,3);
CIf0PertOffsetLP   = CR.audioMf0MeanPertLP(:,4);

meanf0PertOnsetuLP  = CR.audioMf0MeanPertuLP(:,1);
CIf0PertOnsetuLP    = CR.audioMf0MeanPertuLP(:,2);
meanf0PertOffsetuLP = CR.audioMf0MeanPertuLP(:,3);
CIf0PertOffsetuLP   = CR.audioMf0MeanPertuLP(:,4);

meanf0PertOnsetCC  = CR.audioMf0MeanPertCC(:,1);
CIf0PertOnsetCC    = CR.audioMf0MeanPertCC(:,2);
meanf0PertOffsetCC = CR.audioMf0MeanPertCC(:,3);
CIf0PertOffsetCC   = CR.audioMf0MeanPertCC(:,4);

UpR = max([meanf0PertOffsetLP; meanf0PertOffsetuLP; meanf0PertOffsetCC]);
LwR = min([meanf0PertOnsetLP; meanf0PertOnsetuLP; meanf0PertOnsetCC]);

limits = [-0.5 1.0 LwR-30 UpR+30];
% limits = checkLims(limitsM, limitsV);
% pValueThresh = 0.05;

% Plotting Variables
plotpos        = [10 100];
plotdim        = targPixDim;
MeanSubjf0Resp = figure('Color', [1 1 1]);
set(MeanSubjf0Resp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

curSess(strfind(curSess, '_')) = ' ';
dottedStartx = [0 0];
dottedy      = [-500 500];
LPColor      = 'b';
uLPColor     = 'r';
CCColor      = 'g';
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
nM = shadedErrorBar(time, meanf0PertOnsetLP, CIf0PertOnsetLP, 'lineprops', LPColor, 'transparent', 1); %Masked
hold on
nV = shadedErrorBar(time, meanf0PertOnsetuLP, CIf0PertOnsetuLP, 'lineprops', uLPColor, 'transparent', 1); %Voice
hold on
nc = shadedErrorBar(time, meanf0PertOnsetCC, CIf0PertOnsetCC, 'lineprops', CCColor, 'transparent', 1); %Voice

set(nM.mainLine, 'LineWidth', lineThick)
set(nV.mainLine, 'LineWidth', lineThick)
set(nc.mainLine, 'LineWidth', lineThick)
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
fM = shadedErrorBar(time, meanf0PertOffsetLP, CIf0PertOffsetLP, 'lineprops', LPColor, 'transparent', 1); %Masked
hold on
fV = shadedErrorBar(time, meanf0PertOffsetuLP, CIf0PertOffsetuLP, 'lineprops', uLPColor, 'transparent', 1); %Voice
hold on
fc = shadedErrorBar(time, meanf0PertOffsetCC, CIf0PertOffsetCC, 'lineprops', CCColor, 'transparent', 1); %Voice

set(fM.mainLine, 'LineWidth', lineThick)
set(fV.mainLine, 'LineWidth', lineThick)
set(fc.mainLine, 'LineWidth', lineThick)
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
% annoStim = ['SM (M/NM): ' num2str(statSMM) ' cents / ' num2str(statSMV) ' cents'];
% annoResp = ['RM (M/NM): ' num2str(statRMM) ' cents / ' num2str(statRMV) ' cents'];
% annoPerc = ['RP (M/NM): ' num2str(statRPM) '% / ' num2str(statRPV) '%'];
% 
% annoStim = checkSig(statSP, pValueThresh, annoStim);
% annoResp = checkSig(statRP, pValueThresh, annoResp);
% annoPerc = checkSig(statPP, pValueThresh, annoPerc);
%  
% statBox = annotation('textbox',[.30 .75 0.45 0.1],...
%                      'string', {annoStim;
%                                 annoResp
%                                 annoPerc},...
%                       'LineStyle','none',...
%                       'FontName', fontN,...
%                       'FontSize', legAnnoFSize,...
%                       'FontWeight','bold');

legend([fC.mainLine fM.mainLine fV.mainLine fc.mainLine],...
        {[num2str(numCont) ' Control Trials'],...
         [num2str(numPertLP) ' LP Trials'],...
         [num2str(numPertuLP) ' uLP Trials'],...
         [num2str(numPertCC) ' CC Trials']},...
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

function limits = checkLims(limitsM, limitsV)

if limitsM(3) < limitsV(3)
    lwLimit = limitsM(3);
else
    lwLimit = limitsV(3);
end

if limitsM(4) > limitsV(4)
    upLimit = limitsM(4);
else
    upLimit = limitsV(4);
end

limits    = limitsV;
limits(3) = lwLimit;
limits(4) = upLimit;
end

function anno = checkSig(stat, thresh, anno)

if stat < thresh
    anno = [anno '*'];
end
end