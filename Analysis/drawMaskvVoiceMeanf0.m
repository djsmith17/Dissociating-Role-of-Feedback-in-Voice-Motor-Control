function drawMaskvVoiceMeanf0(resM, resV, statLib, targPixDim, pltName, plotFolder)

curSess          = resM.subject;
numCont          = resM.numContTrials;
numMasked        = resM.numPertTrials;
numVoiced        = resV.numPertTrials;

time              = resM.secTime;
meanf0PertOnsetM  = resM.audioMf0MeanPert(:,1);
CIf0PertOnsetM    = resM.audioMf0MeanPert(:,2);
meanf0PertOffsetM = resM.audioMf0MeanPert(:,3);
CIf0PertOffsetM   = resM.audioMf0MeanPert(:,4);

meanf0ContOnsetM  = resM.audioMf0MeanCont(:,1);
CIf0ContOnsetM    = resM.audioMf0MeanCont(:,2);
meanf0ContOffsetM = resM.audioMf0MeanCont(:,3);
CIf0ContOffsetM   = resM.audioMf0MeanCont(:,4);
limitsM           = resM.limitsAmean;

meanf0PertOnsetV  = resV.audioMf0MeanPert(:,1);
CIf0PertOnsetV    = resV.audioMf0MeanPert(:,2);
meanf0PertOffsetV = resV.audioMf0MeanPert(:,3);
CIf0PertOffsetV   = resV.audioMf0MeanPert(:,4);
limitsV           = resV.limitsAmean;

statSMM = round(statLib(1), 1);
statSMV = round(statLib(2), 1);
statRMM = round(statLib(3), 1);
statRMV = round(statLib(4), 1);
statRPM = round(statLib(5));
statRPV = round(statLib(6));
statSP  = statLib(7);
statRP  = statLib(8);
statPP  = statLib(9);

limits = checkLims(limitsM, limitsV);
pValueThresh = 0.05;
figureL      = pltName(end);

% Figure Variables
plotpos        = [10 100];
plotdim        = targPixDim;
IndiSubjf0Resp = figure('Color', [1 1 1]);
set(IndiSubjf0Resp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

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
plot(dottedStartx, dottedy, 'color', [0.3 0.3 0.3], 'LineWidth', lineThick)
hold on
nC = shadedErrorBar(time, meanf0ContOnsetM, CIf0ContOnsetM, 'lineprops', '-k', 'transparent', 1); % Control condition
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
fC = shadedErrorBar(time, meanf0ContOffsetM, CIf0ContOffsetM, 'lineprops', '-k', 'transparent', 1); %Voice
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
            
figureMark = annotation('textbox', [0.01 0.88 0.05 0.1],...
                        'string', figureL,...
                        'LineStyle', 'none',...
                        'FontName', fontN,...
                        'FontSize', titleFSize,...
                        'FontWeight','bold');

legend([fC.mainLine fM.mainLine fV.mainLine],{[num2str(numCont) ' Control Trials'], [num2str(numMasked) ' Masked Trials'], [num2str(numVoiced) ' Not Masked Trials']},...
            'Position', [0.83 0.75 0.1 0.1],...
            'Box', 'off',...
            'Edgecolor', [1 1 1],...
            'FontName', fontN,...
            'FontSize', legAnnoFSize,...
            'FontWeight', 'bold');
        
plots = {'Figure'};
for i = 1:length(plots)
    plTitle = [pltName '.jpg'];
    plTitleT = [pltName '.tif'];

    saveFileName = fullfile(plotFolder, plTitle);
    saveFileNameT= fullfile(plotFolder, plTitleT);
%     X = setImageSaveTif(saveFileNameT, IndiSubjf0Resp);
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

function X = setImageSaveTif(plotFolder, IndiSubjf0Resp)

F = getframe(IndiSubjf0Resp);
[X, Map] = frame2im(F);

imwrite(X, plotFolder, 'Resolution', 300)
end