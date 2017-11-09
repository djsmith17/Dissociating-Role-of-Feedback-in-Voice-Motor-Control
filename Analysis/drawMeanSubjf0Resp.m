function drawMeanSubjf0Resp(allSubjRes, statLib, targPixDim, pltName, plotFolder)

curSess          = 'Mean Participant Response';
numMasked        = allSubjRes.numMaskedTrials;
numVoiced        = allSubjRes.numVoicedTrials;

time              = allSubjRes.secTime;
meanf0PertOnsetM  = allSubjRes.audioMf0MeanPertM(:,1);
CIf0PertOnsetM    = allSubjRes.audioMf0MeanPertM(:,2);
meanf0PertOffsetM = allSubjRes.audioMf0MeanPertM(:,3);
CIf0PertOffsetM   = allSubjRes.audioMf0MeanPertM(:,4);
limitsM           = allSubjRes.limitsAmeanM;

meanf0PertOnsetV  = allSubjRes.audioMf0MeanPertV(:,1);
CIf0PertOnsetV    = allSubjRes.audioMf0MeanPertV(:,2);
meanf0PertOffsetV = allSubjRes.audioMf0MeanPertV(:,3);
CIf0PertOffsetV   = allSubjRes.audioMf0MeanPertV(:,4);
limitsV           = allSubjRes.limitsAmeanV;

statSMM = round(10*statLib(1))/10;
statSMV = round(10*statLib(2))/10;
statRMM = round(10*statLib(3))/10;
statRMV = round(10*statLib(4))/10;
statRPM = round(statLib(5));
statRPV = round(statLib(6));
statSP  = statLib(7);
statRP  = statLib(8);
statPP  = statLib(9);

limits = checkLims(limitsM, limitsV);
pValueThresh = 0.05;

% Plotting Variables
plotpos        = [10 100];
plotdim        = targPixDim;
MeanSubjf0Resp = figure('Color', [1 1 1]);
set(MeanSubjf0Resp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

dottedStartx = [0 0];
dottedy      = [-500 500];
maskColor    = [0 0 1];
voicColor    = [1 0 0];
fontN        = 'Arial';
legAnnoFSize = 30;
titleFSize   = 35;
axisLSize    = 30;
lineThick    = 4;

ha = tight_subplot(1,2,[0.1 0.03],[0.12 0.15],[0.05 0.05]);

%Onset of Perturbation
axes(ha(1))
nM = shadedErrorBar(time, meanf0PertOnsetM, CIf0PertOnsetM, maskColor, 1); %Masked
hold on
nV = shadedErrorBar(time, meanf0PertOnsetV, CIf0PertOnsetV, voicColor, 1); %Voice
hold on
plot(dottedStartx, dottedy,'k','LineWidth',lineThick)

set(nM.mainLine, 'LineWidth', lineThick)
set(nV.mainLine, 'LineWidth', lineThick)
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
fM = shadedErrorBar(time, meanf0PertOffsetM, CIf0PertOffsetM, maskColor, 1); %Masked
hold on
fV = shadedErrorBar(time, meanf0PertOffsetV, CIf0PertOffsetV, voicColor, 1); %Voice
hold on
plot(dottedStartx, dottedy,'k','LineWidth',lineThick)

set(fM.mainLine, 'LineWidth', lineThick)
set(fV.mainLine, 'LineWidth', lineThick)
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
 
statBox = annotation('textbox',[.25 .75 0.45 0.1],...
                     'string', {annoStim;
                                annoResp
                                annoPerc},...
                      'LineStyle','none',...
                      'FontName', fontN,...
                      'FontSize', legAnnoFSize,...
                      'FontWeight','bold');

legend([fM.mainLine fV.mainLine],{[num2str(numMasked) ' Masked Trials'], [num2str(numVoiced) ' Not Masked Trials']},...
            'Position', [0.80 0.75 0.1 0.1],...
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